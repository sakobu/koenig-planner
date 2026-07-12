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
//!     w_meters: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
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
use koenig_damico_planner::solver::{sweep_dual, SweepResult};
use koenig_damico_planner::solver::{sweep_solve as sweep_solve_core, SweepSolveResult};
use koenig_damico_planner::{
    primer_history, solve, solve_from_initial_times, CostModel, PrimerHistory, Solution,
    SolveParams, TimeGrid,
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

/// Monomorphize `solve`/`solve_from_initial_times` over a concrete cost type,
/// then reconstruct the primer-vector history from the converged dual.
///
/// This private helper avoids repeating the dispatch body for each cost
/// variant.  `dyn CostModel` is intentionally not used here: it does not
/// satisfy the `C: CostModel` bound required by the generic `solve` functions.
/// It is also the single point that still holds the concrete `cost` and `grid`
/// (`TimeGrid` is `Copy`, so `solve` consumes a copy), so the primer history is
/// computed here rather than at the response seam.
fn dispatch<C: CostModel>(
    dyn_: &J2Roe,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
    initial_times: Option<&[f64]>,
) -> Result<(Solution, PrimerHistory), PlannerError> {
    let sol = match initial_times {
        Some(ts) => solve_from_initial_times(dyn_, cost, w, grid, params, ts),
        None => solve(dyn_, cost, w, grid, params),
    }?;
    let primer = primer_history(dyn_, cost, &grid, &sol.lambda)?;
    Ok((sol, primer))
}

/// Absolute epoch `[s]` of the chief's first perigee at/after `t_i`, used to
/// place the default eq.-49 perigee windows when the caller omits `t_perigee0`.
///
/// The chief's mean anomaly is anchored at `t_i` (see [`J2Roe`]), and the cost
/// selector compares **absolute** grid times, so the perigee epoch is
/// `t_i + time_to_perigee()`. Omitting `t_i` would shift every FaceMax window by
/// `t_i (mod period)`. Reduces to `time_to_perigee()` for the `t_i = 0` worked
/// example.
///
/// Ref: \[KD20\] eq. 49 (piecewise perigee windows).
fn default_perigee_epoch(chief: &AbsoluteOrbit, t_i: f64) -> f64 {
    t_i + chief.time_to_perigee()
}

/// Resolve the Piecewise cost from a request's optional `period` / `t_perigee0`:
/// default the period to the chief's Keplerian period and the perigee epoch to
/// the p. 10 default (\[KD20\] eq. 49). Shared by [`run`] and [`sweep`].
fn resolve_piecewise(
    chief: &AbsoluteOrbit,
    t_i: f64,
    period: Option<f64>,
    t_perigee0: Option<f64>,
) -> Result<Piecewise, ApiError> {
    let period = period.unwrap_or_else(|| TAU / chief.mean_motion());
    match t_perigee0 {
        Some(tp) => Piecewise::with_perigee_epoch(period, tp),
        None => Piecewise::with_perigee_epoch(period, default_perigee_epoch(chief, t_i)),
    }
    .map_err(bad_request)
}

/// The per-request planning context shared by [`run`] and [`sweep`]: the chief
/// orbit, J2-ROE dynamics, time grid, and resolved solver params. Built once so
/// a sweep over many targets does not reconstruct any of it per target.
struct Context {
    chief: AbsoluteOrbit,
    dyn_: J2Roe,
    grid: TimeGrid,
    params: SolveParams,
}

/// Build the planning context from a request: chief (deg→rad), J2-ROE dynamics,
/// uniform grid, the [`MAX_GRID_POINTS`] guard, and merged params.
fn build_context(req: &SolveRequest) -> Result<Context, ApiError> {
    // 1. Chief mean absolute orbit (angles degrees → radians).
    let chief = AbsoluteOrbit::new(
        req.chief.a,
        req.chief.e,
        req.chief.i.to_radians(),
        req.chief.raan.to_radians(),
        req.chief.argp.to_radians(),
        req.chief.mean_anom.to_radians(),
    );

    // 2. J2-ROE dynamics (validates chief + window).
    let dyn_ = J2Roe::new(chief, req.t_i, req.t_f).map_err(bad_request)?;

    // 3. Uniform time grid (validates dt > 0, t_f > t_i).
    let grid = TimeGrid::uniform(req.t_i, req.t_f, req.dt).map_err(bad_request)?;

    // 3a. Bound the grid size before solving: the Γ-cache allocation and the
    // per-iteration contact sweep are O(grid.len()), driven by attacker-controlled
    // scalars with no upper bound. Reject oversized grids as a bad request before
    // any allocation. The `f64 → usize` saturation in `len()` keeps this correct
    // even for absurd `t_f` (saturates to usize::MAX > cap).
    let n_points = grid.len();
    if n_points > MAX_GRID_POINTS {
        return Err(ApiError {
            kind: ApiErrorKind::BadRequest,
            message: format!(
                "grid has {n_points} points (> {MAX_GRID_POINTS} max); \
                 reduce (t_f - t_i)/dt"
            ),
        });
    }

    // 4. Merge optional parameter overrides with the p. 10-prose defaults.
    let params = resolve_params(req.params.clone());

    Ok(Context {
        chief,
        dyn_,
        grid,
        params,
    })
}

// ──────────────────────────────────────────────────────────────────────────────
// Public entry point
// ──────────────────────────────────────────────────────────────────────────────

/// Maximum accepted request-body size in bytes for [`run_json`]. The uncapped
/// library/py/wasm entrypoints have no transport-layer limit, so this bounds
/// the worst-case parse allocation. The HTTP server applies its own 64 KiB cap.
pub const MAX_REQUEST_BYTES: usize = 1_048_576;

/// Maximum number of grid points [`run`] will solve over.
///
/// The largest real request (the worked example) is ~3,934 points; this is ~25×
/// that, bounding the Γ-cache (`grid.len() × 144 B`) to ~14 MB and a solve to tens
/// of ms while never rejecting a realistic mission horizon. This is the
/// untrusted-boundary guard against the grid-size complexity DoS: all
/// three frontends (HTTP, WASM, Python) funnel through [`run`], so this single cap
/// protects every one.
pub const MAX_GRID_POINTS: usize = 100_000;

/// Maximum number of targets [`sweep`] / [`sweep_solve`] evaluate in one batch.
///
/// Both run `w_list.len()` solves — [`sweep`] a dual solve each, [`sweep_solve`]
/// a full primal [`solve`] each — over the (already
/// [`MAX_GRID_POINTS`]-bounded) window, so `w_list` is the second
/// attacker-controlled cost dimension. This bounds it with the same discipline
/// as the grid: a reachable-set trace needs ~180 targets and a dense cost-map
/// grid a few thousand, so this never rejects a realistic batch.
pub const MAX_SWEEP_TARGETS: usize = 100_000;

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
    let Context {
        chief,
        dyn_,
        grid,
        params,
    } = build_context(&req)?;

    // Nondimensionalize the target pseudostate (divide by chief.a).
    let w = Pseudostate::from_row_slice(&req.w_meters) / chief.a;

    // 6. Dispatch per cost model (monomorphize per match arm). Each arm returns
    //    the solution paired with its primer-vector history.
    let its = req.initial_times.as_deref();
    let (sol, primer) = match req.cost {
        CostSpec::Norm2 => dispatch(&dyn_, &ConstNorm2(Norm2), w, grid, &params, its),
        CostSpec::FaceMax => dispatch(&dyn_, &ConstFaceMax(FaceMax), w, grid, &params, its),
        CostSpec::Piecewise { period, t_perigee0 } => {
            let cost = resolve_piecewise(&chief, req.t_i, period, t_perigee0)?;
            dispatch(&dyn_, &cost, w, grid, &params, its)
        }
    }
    .map_err(map_dispatch_error)?;

    // 7. Finite-guard: serde_json renders non-finite f64 as `null`. Covers the
    //    solution and the primer-vector history (magnitudes and RTN components).
    if !sol.total_dv.is_finite()
        || !sol.residual.is_finite()
        || sol
            .maneuvers
            .iter()
            .any(|m| !m.dv.iter().all(|x| x.is_finite()))
        || !sol.lambda.iter().all(|x| x.is_finite())
        || !primer.magnitudes.iter().all(|g| g.is_finite())
        || primer
            .vectors
            .iter()
            .any(|p| !p.iter().all(|x| x.is_finite()))
    {
        return Err(ApiError {
            kind: ApiErrorKind::Solver,
            message: "solver produced a non-finite result".into(),
        });
    }

    // 8. Map (Solution, PrimerHistory) → SolveResponse via the field-exhaustive
    //    `From` in convert.rs: a new field on either becomes a compile error here.
    Ok((sol, primer).into())
}

/// Parse a JSON [`SolveRequest`], run it, and serialize the [`SolveResponse`] to JSON.
///
/// The shared string-in / string-out entry point reused by the WASM and Python
/// frontends so the serde glue lives in exactly one place.
///
/// # Errors
/// Returns [`ApiError`] with `kind = "bad_request"` for malformed request JSON,
/// invalid inputs, or a body exceeding [`MAX_REQUEST_BYTES`], or `kind = "solver"`
/// for numerically unsolvable / failed problems (including an internal
/// response-serialization failure).
pub fn run_json(input: &str) -> Result<String, ApiError> {
    if input.len() > MAX_REQUEST_BYTES {
        return Err(ApiError {
            kind: ApiErrorKind::BadRequest,
            message: format!(
                "request body is {} bytes (> {MAX_REQUEST_BYTES} max)",
                input.len()
            ),
        });
    }
    let req: SolveRequest = serde_json::from_str(input).map_err(|e| ApiError {
        kind: ApiErrorKind::BadRequest,
        message: format!("invalid request JSON: {e}"),
    })?;
    let resp = run(req)?;
    serde_json::to_string(&resp).map_err(|e| ApiError {
        kind: ApiErrorKind::Solver,
        message: format!("failed to serialize response: {e}"),
    })
}

/// Evaluate the min-fuel dual gauge for many targets over `base`'s window.
///
/// Builds the chief/dynamics/grid/cost once from `base`, nondimensionalizes each
/// `w_list` entry (meters ÷ `chief.a`), and returns one [`SweepPoint`] per
/// target: the gauge `c*` (m/s; `None` if unreachable) and the dual normal `λ`.
/// This is the batch sibling of [`run`]; it never returns maneuvers.
///
/// # Errors
/// Returns [`ApiError`] with `kind = "bad_request"` for an invalid orbit, a
/// degenerate or oversized grid ([`MAX_GRID_POINTS`]), or more than
/// [`MAX_SWEEP_TARGETS`] targets. Per-target unreachability is reported as
/// `SweepPoint { feasible: false, c_star: None, .. }`, not an error.
pub fn sweep(base: &SolveRequest, w_list: &[[f64; 6]]) -> Result<Vec<SweepPoint>, ApiError> {
    let ctx = build_context(base)?;

    // Bound the batch size: `sweep` runs one dual solve per target, an
    // attacker-controlled count orthogonal to the grid guard in `build_context`.
    if w_list.len() > MAX_SWEEP_TARGETS {
        return Err(ApiError {
            kind: ApiErrorKind::BadRequest,
            message: format!(
                "sweep has {} targets (> {MAX_SWEEP_TARGETS} max); reduce w_list",
                w_list.len()
            ),
        });
    }

    let targets: Vec<Pseudostate> = w_list
        .iter()
        .map(|w| Pseudostate::from_row_slice(w) / ctx.chief.a)
        .collect();

    // Monomorphize sweep_dual per cost model, exactly as `run`/`dispatch` do.
    let results = match base.cost {
        CostSpec::Norm2 => sweep_dual(&ctx.dyn_, &ConstNorm2(Norm2), &ctx.grid, &targets),
        CostSpec::FaceMax => sweep_dual(&ctx.dyn_, &ConstFaceMax(FaceMax), &ctx.grid, &targets),
        CostSpec::Piecewise { period, t_perigee0 } => {
            let cost = resolve_piecewise(&ctx.chief, base.t_i, period, t_perigee0)?;
            sweep_dual(&ctx.dyn_, &cost, &ctx.grid, &targets)
        }
    }
    .map_err(map_dispatch_error)?;

    Ok(results.into_iter().map(sweep_point).collect())
}

/// Map a core [`SweepResult`] to a wire [`SweepPoint`], scrubbing non-finite /
/// infeasible results to `c_star: None` and `lambda: [0; 6]` (serde renders a
/// non-finite f64 as `null`, a type lie inside a number array).
fn sweep_point(r: SweepResult) -> SweepPoint {
    // Destructure with no `..` so a new SweepResult field must be handled here.
    let SweepResult {
        c_star,
        lambda,
        feasible,
    } = r;
    let finite = feasible && c_star.is_finite() && lambda.iter().all(|x| x.is_finite());
    let mut lambda_out = [0.0_f64; 6];
    if finite {
        lambda_out.copy_from_slice(lambda.as_slice());
    }
    SweepPoint {
        c_star: finite.then_some(c_star),
        lambda: lambda_out,
        feasible: finite,
    }
}

/// Evaluate the min-fuel primal for many targets over `base`'s window.
///
/// Builds the chief/dynamics/grid/params once from `base` (as [`sweep`] does),
/// nondimensionalizes each `w_list` entry, and returns one [`SweepSolvePoint`]
/// per target: cost, dual, Algorithm 2 confidence, and burn count. Primal
/// batch sibling of [`run`]; unlike [`sweep`] (the exact dual gauge via
/// [`crate::core::solver::refine_socp()`]), this loops the full active-set
/// [`solve`] per target, so it never returns maneuvers but does return a
/// confidence signal.
///
/// # Errors
/// Returns [`ApiError`] with `kind = "bad_request"` for an invalid orbit, a
/// degenerate or oversized grid ([`MAX_GRID_POINTS`]), or more than
/// [`MAX_SWEEP_TARGETS`] targets. Per-target unreachability is reported as
/// `SweepSolvePoint { feasible: false, c_star: None, .. }`, not an error.
pub fn sweep_solve(base: &SolveRequest, w_list: &[[f64; 6]]) -> Result<Vec<SweepSolvePoint>, ApiError> {
    let ctx = build_context(base)?;

    // Bound the batch size, same discipline as `sweep`.
    if w_list.len() > MAX_SWEEP_TARGETS {
        return Err(ApiError {
            kind: ApiErrorKind::BadRequest,
            message: format!(
                "sweep_solve has {} targets (> {MAX_SWEEP_TARGETS} max); reduce w_list",
                w_list.len()
            ),
        });
    }

    let targets: Vec<Pseudostate> = w_list
        .iter()
        .map(|w| Pseudostate::from_row_slice(w) / ctx.chief.a)
        .collect();

    // `build_context` already resolved params into `ctx.params` (via
    // `resolve_params`), exactly as `run` consumes them — use it directly.
    let results = match base.cost {
        CostSpec::Norm2 => sweep_solve_core(&ctx.dyn_, &ConstNorm2(Norm2), &ctx.grid, &targets, &ctx.params),
        CostSpec::FaceMax => sweep_solve_core(&ctx.dyn_, &ConstFaceMax(FaceMax), &ctx.grid, &targets, &ctx.params),
        CostSpec::Piecewise { period, t_perigee0 } => {
            let cost = resolve_piecewise(&ctx.chief, base.t_i, period, t_perigee0)?;
            sweep_solve_core(&ctx.dyn_, &cost, &ctx.grid, &targets, &ctx.params)
        }
    };

    Ok(results.into_iter().map(sweep_solve_point).collect())
}

/// Map a core [`SweepSolveResult`] to a wire [`SweepSolvePoint`], scrubbing
/// non-finite / infeasible results to `c_star: None`, `residual: None`,
/// `lambda: [0; 6]`, and `n_maneuvers: 0` — mirrors [`sweep_point`].
fn sweep_solve_point(r: SweepSolveResult) -> SweepSolvePoint {
    // Destructure with no `..` so a new SweepSolveResult field must be handled here.
    let SweepSolveResult {
        c_star,
        lambda,
        feasible,
        iterations,
        residual,
        n_maneuvers,
    } = r;
    let finite = feasible && c_star.is_finite() && lambda.iter().all(|x| x.is_finite());
    let mut lambda_out = [0.0_f64; 6];
    if finite {
        lambda_out.copy_from_slice(lambda.as_slice());
    }
    SweepSolvePoint {
        c_star: finite.then_some(c_star),
        lambda: lambda_out,
        feasible: finite,
        iterations,
        residual: (finite && residual.is_finite()).then_some(residual),
        n_maneuvers: if finite { n_maneuvers } else { 0 },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use koenig_damico_planner::cost::Piecewise;

    // Ref: [KD20] eq. 49 — the default perigee windows must sit on the chief's
    // ACTUAL perigee. The chief's mean anomaly is anchored at t_i (J2Roe) and the
    // cost selector compares absolute grid times, so the absolute perigee epoch is
    // t_i + time_to_perigee(). A default omitting t_i misplaces the FaceMax window
    // by t_i (mod period) for any t_i != 0.
    #[test]
    fn default_perigee_epoch_places_window_on_true_perigee_for_nonzero_t_i() {
        let chief = AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            90.0_f64.to_radians(),
        );
        let period = TAU / chief.mean_motion();
        // A t_i whose remainder mod period exceeds the 3600 s window half-width,
        // so omitting it moves the true perigee clear out of the window.
        let t_i = 10_000.0;
        let true_perigee = t_i + chief.time_to_perigee(); // absolute time, M ≡ 0

        assert_eq!(default_perigee_epoch(&chief, t_i), true_perigee);

        // FaceMax is active at the true perigee under the correct (t_i-anchored)
        // epoch...
        let correct = Piecewise::with_perigee_epoch(period, true_perigee).unwrap();
        assert!(correct.in_perigee_window(true_perigee));
        // ...but the t_i-less epoch (the pre-fix default) puts the window a full
        // t_i off, so the true perigee falls outside it.
        let buggy = Piecewise::with_perigee_epoch(period, chief.time_to_perigee()).unwrap();
        assert!(!buggy.in_perigee_window(true_perigee));
    }

    fn worked_example_request() -> SolveRequest {
        SolveRequest {
            chief: OrbitDto {
                a: 25_000e3,
                e: 0.7,
                i: 40.0,
                raan: 358.0,
                argp: 0.0,
                mean_anom: 180.0,
            },
            t_i: 0.0,
            t_f: 117_990.0,
            dt: 300.0, // coarse grid keeps the test fast; still richly reachable
            w_meters: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
            cost: CostSpec::Piecewise {
                period: None,
                t_perigee0: None,
            },
            params: None,
            initial_times: None,
        }
    }

    // Weak duality guarantees the dual gauge c* ≤ run()'s primal total_dv on the
    // same grid. Strong duality (SOCP + Slater) makes the true gap ~solver
    // tolerance; the few-% slack below absorbs run()'s active-set / extract
    // suboptimality. Both targets are reachable.
    #[test]
    fn sweep_matches_run_within_duality_gap() {
        let base = worked_example_request();
        let w_a = base.w_meters;
        let w_b = [25.0, 2500.0, 50.0, 50.0, 0.0, 200.0];

        let pts = sweep(&base, &[w_a, w_b]).unwrap();
        assert_eq!(pts.len(), 2);
        assert!(pts[0].feasible && pts[1].feasible);

        let primal_a = run(base.clone()).unwrap().total_dv;
        let dual_a = pts[0].c_star.unwrap();
        assert!(
            dual_a <= primal_a + 1e-9,
            "dual {dual_a} should be ≤ primal {primal_a}"
        );
        assert!(
            (primal_a - dual_a) <= 0.05 * primal_a,
            "duality gap too large: dual {dual_a} vs primal {primal_a}"
        );
    }

    // Infeasible / non-finite results scrub to c_star: None and lambda: [0; 6].
    #[test]
    fn sweep_point_scrubs_infeasible_to_none() {
        let infeasible = SweepResult {
            c_star: f64::NAN,
            lambda: nalgebra::SVector::<f64, 6>::zeros(),
            feasible: false,
        };
        assert_eq!(sweep_point(infeasible).c_star, None);

        // A feasible flag but a non-finite lambda component: still scrubbed, and
        // the non-finite lambda is zeroed rather than leaked as a JSON null.
        let mut bad_lambda = nalgebra::SVector::<f64, 6>::zeros();
        bad_lambda[0] = f64::INFINITY;
        let nonfinite = SweepResult {
            c_star: 1.0,
            lambda: bad_lambda,
            feasible: true,
        };
        let p = sweep_point(nonfinite);
        assert_eq!(p.c_star, None);
        assert!(!p.feasible);
        assert_eq!(p.lambda, [0.0; 6]);
    }

    // w_list beyond MAX_SWEEP_TARGETS is a bad request, not a huge batch: the
    // guard rejects before any solve (mirrors the MAX_GRID_POINTS grid guard).
    #[test]
    fn sweep_rejects_oversized_w_list() {
        let base = worked_example_request();
        let too_many = vec![[0.0; 6]; MAX_SWEEP_TARGETS + 1];
        let err = sweep(&base, &too_many).unwrap_err();
        assert_eq!(err.kind, ApiErrorKind::BadRequest);
    }

    // Single-grid-time window (t_f - t_i < dt ⇒ TimeGrid::uniform yields exactly
    // one grid time) so a generic 6-D target is structurally unreachable, same
    // fixture shape as the core sweep_solve's unreachable_target_is_infeasible.
    fn short_window_request() -> SolveRequest {
        SolveRequest {
            t_i: 0.0,
            t_f: 10.0,
            dt: 100.0,
            ..worked_example_request()
        }
    }

    // sweep_solve calls the identical solve() with the identical resolved
    // params as run(), so its c_star must equal run()'s total_dv exactly (not
    // just within the sweep_dual duality gap).
    #[test]
    fn sweep_solve_matches_run() {
        let base = worked_example_request();
        let w = base.w_meters;
        let run_dv = run(base.clone()).unwrap().total_dv;
        let pts = sweep_solve(&base, &[w]).unwrap();
        assert_eq!(pts.len(), 1);
        assert!(pts[0].feasible);
        let c = pts[0].c_star.unwrap();
        assert!((c - run_dv).abs() <= 1e-9, "sweep_solve c*={c} vs run={run_dv}");
        assert!(pts[0].n_maneuvers >= 1);
    }

    // Single-grid-time window → generic 6-D target structurally unreachable →
    // solve()'s internal dual is unbounded → Err → scrubbed to feasible:false.
    #[test]
    fn sweep_solve_scrubs_infeasible_to_none() {
        let base = short_window_request();
        let pts = sweep_solve(&base, &[[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]).unwrap();
        assert!(!pts[0].feasible);
        assert!(pts[0].c_star.is_none());
        assert!(pts[0].residual.is_none());
    }

    // w_list beyond MAX_SWEEP_TARGETS is a bad request, mirroring sweep_rejects_oversized_w_list.
    #[test]
    fn sweep_solve_rejects_oversized_w_list() {
        let base = worked_example_request();
        let too_many = vec![[0.0; 6]; MAX_SWEEP_TARGETS + 1];
        assert!(sweep_solve(&base, &too_many).is_err());
    }

    /// The worked example's near-equatorial twin: identical HEO chief, window,
    /// cost, and target, but the inclination is pulled from 40° to 0.001°
    /// (~1.7e-5 rad — well clear of the 1e-9 rad singularity guard). Per
    /// `J2Roe`'s `# Conditioning` note, that near-equatorial regime blows the
    /// cross-track `B(t)` entries up to ~1e5 (vs ~1e-4 typical), so the reachable
    /// set is extremely anisotropic — ill-conditioned — while the target stays
    /// COMFORTABLY reachable (feasible, finite c*, several maneuvers), not a
    /// reachability-boundary probe.
    fn near_equatorial_request() -> SolveRequest {
        let base = worked_example_request();
        SolveRequest {
            chief: OrbitDto {
                i: 0.001,
                ..base.chief
            },
            ..base
        }
    }

    // Validates the frontend's per-cell "confidence" fields on an ill-conditioned
    // BUT comfortably-reachable cell (both cells feasible ⇒ not a boundary probe).
    //
    // FINDING: `residual` is the signal that discriminates conditioning — the
    // relative recovery error jumps ~7 orders of magnitude, from ~9e-15 on the
    // clean i=40° cell to ~2e-7 on the near-equatorial twin. The refine
    // `iterations` count does NOT track conditioning: it wanders in a small
    // [1,5] band and here is actually LOWER on the ill-conditioned cell (2 vs 3),
    // so the `moved` OR below reduces to the residual term — i.e. this doubles as
    // a residual-regression guard (if residual collapsed back to ~machine-zero,
    // the iterations term could not rescue the assertion).
    #[test]
    fn confidence_signal_moves_on_an_ill_conditioned_target() {
        let clean_req = worked_example_request();
        let clean = sweep_solve(&clean_req, &[clean_req.w_meters])
            .unwrap()
            .remove(0);
        let marg_req = near_equatorial_request();
        let marginal = sweep_solve(&marg_req, &[marg_req.w_meters])
            .unwrap()
            .remove(0);

        // Both reachable — this is not a reachability-boundary probe.
        assert!(clean.feasible && marginal.feasible);

        let moved = marginal.iterations > clean.iterations
            || marginal.residual.unwrap_or(0.0) > clean.residual.unwrap_or(0.0) * 5.0;
        assert!(
            moved,
            "confidence signal did not discriminate: clean(iters={}, res={:?}) vs marginal(iters={}, res={:?})",
            clean.iterations, clean.residual, marginal.iterations, marginal.residual
        );
    }
}
