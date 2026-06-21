//! # koenig-damico-planner-api
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
//! use koenig_damico_planner_api::{run, CostSpec, OrbitDto, SolveRequest};
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

mod convert;
mod dto;

pub use dto::*;

use convert::{bad_request, map_dispatch_error, resolve_params};
use koenig_damico_planner::cost::{FaceMax, Norm2, Piecewise, SublevelSet};
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::{
    solve, solve_from_initial_times, CostModel, Solution, SolveParams, TimeGrid,
};
use koenig_damico_planner::{PlannerError, Pseudostate};
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

/// Maximum number of grid points [`run`] will solve over.
///
/// The largest real request (the worked example) is ~3,934 points; this is ~25×
/// that, bounding the Γ-cache (`grid.len() × 144 B`) to ~14 MB and a solve to tens
/// of ms while never rejecting a realistic mission horizon. This is the
/// untrusted-boundary guard against the grid-size complexity DoS (audit B1): all
/// three frontends (HTTP, WASM, Python) funnel through [`run`], so this single cap
/// protects every one.
pub const MAX_GRID_POINTS: usize = 100_000;

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

    // 3a. Bound the grid size before solving (audit B1): the Γ-cache allocation
    // and the per-iteration contact sweep are O(grid.len()), driven by
    // attacker-controlled scalars with no upper bound. Reject oversized grids as a
    // bad request *before* any allocation. The `f64 → usize` saturation in `len()`
    // keeps this correct even for absurd `t_f` (saturates to usize::MAX > cap).
    let n_points = grid.len();
    if n_points > MAX_GRID_POINTS {
        return Err(ApiError {
            kind: "bad_request",
            message: format!(
                "grid has {n_points} points (> {MAX_GRID_POINTS} max); \
                 reduce (t_f - t_i)/dt"
            ),
        });
    }

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

/// Parse a JSON [`SolveRequest`], run it, and serialize the [`SolveResponse`] to JSON.
///
/// The shared string-in / string-out entry point reused by the WASM and Python
/// frontends so the serde glue lives in exactly one place.
///
/// # Errors
/// Returns [`ApiError`] with `kind = "bad_request"` for malformed request JSON
/// or invalid inputs, or `kind = "solver"` for numerically unsolvable / failed
/// problems (including an internal response-serialization failure).
pub fn run_json(input: &str) -> Result<String, ApiError> {
    let req: SolveRequest = serde_json::from_str(input).map_err(|e| ApiError {
        kind: "bad_request",
        message: format!("invalid request JSON: {e}"),
    })?;
    let resp = run(req)?;
    serde_json::to_string(&resp).map_err(|e| ApiError {
        kind: "solver",
        message: format!("failed to serialize response: {e}"),
    })
}
