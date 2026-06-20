//! # koenig-planner-api
//!
//! Serde request/response DTOs and the single [`run`] dispatch entry point for
//! the Koenig-D'Amico maneuver planner.
//!
//! `run` is the **one** place the generic `solve`/`solve_from_initial_times`
//! is monomorphized over the cost model.  HTTP, WASM, and Python frontends
//! all call it and never need to know which cost model was selected.
//!
//! ## Usage
//! ```no_run
//! use koenig_planner_api::{run, CostSpec, OrbitDto, SolveRequest};
//!
//! let req = SolveRequest {
//!     chief: OrbitDto { a: 25_000e3, e: 0.7, i: 40.0, raan: 358.0, argp: 0.0, mean_anom: 180.0 },
//!     t_i: 0.0,
//!     t_f: 117_990.0,
//!     dt: 30.0,
//!     w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
//!     cost: CostSpec::Piecewise { period: None, t_perigee0: None },
//!     params: None,
//!     initial_times: None,
//! };
//! let resp = run(req).expect("worked example should solve");
//! assert!(resp.total_dv > 0.0);
//! ```

// Re-export the core crate so downstream crates do not pin it independently.
pub use koenig_damico_planner as core;

use koenig_damico_planner::cost::{FaceMax, Norm2, Piecewise, SublevelSet};
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::{
    solve, solve_from_initial_times, CostModel, Solution, SolveParams, TimeGrid,
};
use koenig_damico_planner::{PlannerError, Pseudostate};
use serde::{Deserialize, Serialize};
use std::f64::consts::TAU;

// ──────────────────────────────────────────────────────────────────────────────
// Private constant-cost adapters
//
// `Norm2` and `FaceMax` implement `SublevelSet`, not `CostModel`.  The
// algorithm requires `C: CostModel`, so we wrap each in a trivial adapter that
// returns the same sublevel set for every time.  These are private; callers
// select via `CostSpec`.
// ──────────────────────────────────────────────────────────────────────────────

/// Constant `Norm2` cost model: returns `Norm2` for every time.
struct ConstNorm2(Norm2);

impl CostModel for ConstNorm2 {
    fn at(&self, _t: f64) -> &dyn SublevelSet {
        &self.0
    }
}

/// Constant `FaceMax` cost model: returns `FaceMax` for every time.
struct ConstFaceMax(FaceMax);

impl CostModel for ConstFaceMax {
    fn at(&self, _t: f64) -> &dyn SublevelSet {
        &self.0
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Request DTOs
// ──────────────────────────────────────────────────────────────────────────────

/// Chief orbit definition.  Angles are in **degrees** (converted to radians
/// server-side); `a` is in **metres**.
///
/// These are the six mean Keplerian elements `[a, e, i, Ω, ω, M]` as used
/// throughout Koenig & D'Amico (2020).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitDto {
    /// Semimajor axis `[m]`.
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination `[deg]`.
    pub i: f64,
    /// Right ascension of the ascending node, Ω `[deg]`.
    pub raan: f64,
    /// Argument of perigee, ω `[deg]`.
    pub argp: f64,
    /// Mean anomaly, M `[deg]`.
    pub mean_anom: f64,
}

/// Which cost model to apply at each maneuver time.
///
/// For `Piecewise`, `period` defaults to the chief's Keplerian orbit period
/// (`2π / n`) when omitted — supplying a period unrelated to the chief
/// silently misaligns the perigee windows, so prefer the default.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum CostSpec {
    /// L2 norm (isotropic cost).
    Norm2,
    /// FaceMax gauge (fuel-optimal for an impulsive thruster set).
    FaceMax,
    /// Piecewise eq.-49 selector: FaceMax near perigee, Norm2 elsewhere.
    Piecewise {
        /// Orbit period `[s]`.  When `None`, derived as `2π / n` from the
        /// chief — strongly preferred so the perigee windows align correctly.
        #[serde(default)]
        period: Option<f64>,
        /// First perigee-passage epoch `[s]`.  When `None`, defaults to
        /// `period / 2` (apogee-at-`t = 0` assumption, matching the worked
        /// example where `M₀ = 180°`).
        #[serde(default)]
        t_perigee0: Option<f64>,
    },
}

/// Solver tuning knobs.  Every field is optional; missing fields fall back to
/// [`SolveParams::default`] (Table III of Koenig & D'Amico 2020).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolveParamsDto {
    /// Coarse-sample count `|T^d|` for Algorithm 1 (default 20).
    #[serde(default)]
    pub n_coarse: Option<usize>,
    /// Initial candidate-time count `n_init` (default 6).
    #[serde(default)]
    pub n_init: Option<usize>,
    /// Convergence tolerance `eps_cost` (default 0.01).
    #[serde(default)]
    pub eps_cost: Option<f64>,
    /// Slack-removal tolerance `eps_remove` (default 0.01).
    #[serde(default)]
    pub eps_remove: Option<f64>,
}

/// A full planning request.
///
/// Angles in [`OrbitDto`] are **degrees**; `w_metres` is in **metres**;
/// times are in **seconds**.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolveRequest {
    /// Chief mean absolute orbit (angles in degrees).
    pub chief: OrbitDto,
    /// Initial time `t_i` `[s]`.
    pub t_i: f64,
    /// Final time `t_f` `[s]`.
    pub t_f: f64,
    /// Grid step `dt` `[s]`.
    pub dt: f64,
    /// Target pseudostate in **metres**.  The server divides each component
    /// by `chief.a` to produce the dimensionless `w` passed to the planner,
    /// matching the nondimensionalisation in the worked example.
    pub w_metres: [f64; 6],
    /// Cost model selection.
    pub cost: CostSpec,
    /// Optional solver tuning (default = Table III).
    #[serde(default)]
    pub params: Option<SolveParamsDto>,
    /// Optional explicit initial candidate times for Algorithm 2 (bypasses
    /// Algorithm 1 when provided, enabling the paper's initialization study).
    #[serde(default)]
    pub initial_times: Option<Vec<f64>>,
}

// ──────────────────────────────────────────────────────────────────────────────
// Response DTOs
// ──────────────────────────────────────────────────────────────────────────────

/// A single impulsive maneuver in the RTN frame.
#[derive(Debug, Clone, Serialize)]
pub struct ManeuverDto {
    /// Application time `[s]`, measured from `t_i`.
    pub t: f64,
    /// Delta-v `[m/s]`, RTN components `[R, T, N]`.
    pub dv: [f64; 3],
}

/// Successful planner output.
#[derive(Debug, Clone, Serialize)]
pub struct SolveResponse {
    /// Ordered list of maneuvers.
    pub maneuvers: Vec<ManeuverDto>,
    /// Total fuel cost `Σ ‖Δvⱼ‖` `[m/s]`.
    pub total_dv: f64,
    /// Algorithm 2 iteration count.
    pub iterations: usize,
    /// Relative residual `‖w_err‖ / ‖w‖` of the pre-prune min-fuel solution.
    pub residual: f64,
    /// Optimal dual variable `λ_opt ∈ ℝ⁶`.
    pub lambda: [f64; 6],
}

/// Owned error that decouples the wire contract from [`PlannerError`].
///
/// `kind` is the status class for HTTP frontends:
/// - `"bad_request"` — invalid input (the caller should fix the request).
/// - `"solver"` — well-formed input but numerically unsolvable / failed.
#[derive(Debug, thiserror::Error)]
#[error("{kind}: {message}")]
pub struct ApiError {
    /// Status class: `"bad_request"` or `"solver"`.
    pub kind: &'static str,
    /// Human-readable description of what went wrong.
    pub message: String,
}

// ──────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ──────────────────────────────────────────────────────────────────────────────

/// Map a [`PlannerError`] to a `"bad_request"` [`ApiError`].
fn bad_request(e: PlannerError) -> ApiError {
    ApiError {
        kind: "bad_request",
        message: e.to_string(),
    }
}

/// Map a [`PlannerError`] from the dispatch (solve/solve_from_initial_times)
/// to the correct [`ApiError`] kind.
///
/// `InvalidInput` is caller-fixable → `"bad_request"`.
/// All other variants indicate a numerically unsolvable / failed problem →
/// `"solver"`.  The match is exhaustive so a future new variant forces an
/// explicit decision here at compile time.
fn map_dispatch_error(e: koenig_damico_planner::PlannerError) -> ApiError {
    use koenig_damico_planner::PlannerError;
    match e {
        // Caller-fixable: bad request.
        PlannerError::InvalidInput(_) => ApiError {
            kind: "bad_request",
            message: e.to_string(),
        },
        // Well-formed request, numerically unsolvable / solver failure.
        PlannerError::SolverFailed(_)
        | PlannerError::NotConverged { .. }
        | PlannerError::KeplerDivergence { .. } => ApiError {
            kind: "solver",
            message: e.to_string(),
        },
    }
}

/// Merge optional overrides from [`SolveParamsDto`] into the Table III defaults.
fn resolve_params(dto: Option<SolveParamsDto>) -> SolveParams {
    let mut p = SolveParams::default();
    if let Some(d) = dto {
        if let Some(v) = d.n_coarse {
            p.n_coarse = v;
        }
        if let Some(v) = d.n_init {
            p.n_init = v;
        }
        if let Some(v) = d.eps_cost {
            p.eps_cost = v;
        }
        if let Some(v) = d.eps_remove {
            p.eps_remove = v;
        }
    }
    p
}

/// Monomorphize `solve`/`solve_from_initial_times` over a concrete cost type.
///
/// This private helper avoids repeating the dispatch body for each cost
/// variant.  `dyn CostModel` is intentionally not used here: it does not
/// satisfy the `C: CostModel` bound required by the generic `solve` functions.
fn dispatch<C: CostModel>(
    dyn_: &J2Roe,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
    initial_times: Option<&[f64]>,
) -> Result<Solution, PlannerError> {
    match initial_times {
        Some(ts) => solve_from_initial_times(dyn_, cost, w, grid, params, ts),
        None => solve(dyn_, cost, w, grid, params),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Public entry point
// ──────────────────────────────────────────────────────────────────────────────

/// Plan a maneuver set from a serde request.
///
/// This is the **one** place `solve`/`solve_from_initial_times` is
/// monomorphized over the cost model; the HTTP, WASM, and Python frontends
/// all call this function and never touch the generic core API directly.
///
/// # Errors
/// Returns [`ApiError`] with `kind = "bad_request"` for invalid inputs (bad
/// orbit, degenerate grid, …) or `kind = "solver"` for numerically unsolvable
/// / failed problems.
pub fn run(req: SolveRequest) -> Result<SolveResponse, ApiError> {
    // 1. Build the chief mean absolute orbit (angles degrees → radians).
    let chief = AbsoluteOrbit::new(
        req.chief.a,
        req.chief.e,
        req.chief.i.to_radians(),
        req.chief.raan.to_radians(),
        req.chief.argp.to_radians(),
        req.chief.mean_anom.to_radians(),
    );

    // 2. Build the J2-ROE dynamics (validates chief + window).
    let dyn_ = J2Roe::new(chief, req.t_i, req.t_f).map_err(bad_request)?;

    // 3. Build the uniform time grid (validates dt > 0, t_f >= t_i).
    let grid = TimeGrid::uniform(req.t_i, req.t_f, req.dt).map_err(bad_request)?;

    // 4. Nondimensionalize the target pseudostate (divide by chief.a).
    let w = Pseudostate::from_row_slice(&req.w_metres) / chief.a;

    // 5. Merge optional parameter overrides with Table III defaults.
    let params = resolve_params(req.params);

    // 6. Dispatch per cost model (monomorphize per match arm).
    let its = req.initial_times.as_deref();
    let sol = match req.cost {
        CostSpec::Norm2 => dispatch(&dyn_, &ConstNorm2(Norm2), w, grid, &params, its),
        CostSpec::FaceMax => dispatch(&dyn_, &ConstFaceMax(FaceMax), w, grid, &params, its),
        CostSpec::Piecewise { period, t_perigee0 } => {
            let period = period.unwrap_or_else(|| TAU / chief.mean_motion());
            let cost = match t_perigee0 {
                Some(tp) => Piecewise::with_perigee_epoch(period, tp),
                None => Piecewise::new(period),
            };
            dispatch(&dyn_, &cost, w, grid, &params, its)
        }
    }
    .map_err(map_dispatch_error)?;

    // 7. Finite-guard: serde_json renders non-finite f64 as `null`.
    if !sol.total_dv.is_finite()
        || !sol.residual.is_finite()
        || sol
            .maneuvers
            .iter()
            .any(|m| !m.dv.iter().all(|x| x.is_finite()))
        || !sol.lambda.iter().all(|x| x.is_finite())
    {
        return Err(ApiError {
            kind: "solver",
            message: "solver produced a non-finite result".into(),
        });
    }

    // 8. Map Solution → SolveResponse.
    let maneuvers = sol
        .maneuvers
        .iter()
        .map(|m| ManeuverDto {
            t: m.t,
            dv: [m.dv[0], m.dv[1], m.dv[2]],
        })
        .collect();

    let lam = sol.lambda;
    let lambda = [lam[0], lam[1], lam[2], lam[3], lam[4], lam[5]];

    Ok(SolveResponse {
        maneuvers,
        total_dv: sol.total_dv,
        iterations: sol.iterations,
        residual: sol.residual,
        lambda,
    })
}
