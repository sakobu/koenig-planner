//! Orchestration of the three-step algorithm (Init -> Refine -> Extract).

mod extract;
mod init;
mod refine;

use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{Dual, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N};
use nalgebra::SMatrix;

/// Contact value `g_{U(1,t)}(Γᵀ(t)·lambda)` at grid index `k`.
///
/// `gammas[k] = Γ(grid.time(k))` is the 6×3 dynamics matrix; the cost methods
/// operate on the control-space projection `Γᵀ(t)·lambda ∈ ℝ³`.
fn contact_at<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    k: usize,
    lambda: &Dual,
) -> f64 {
    let y = gammas[k].transpose() * lambda;
    cost.at(grid.time(k)).contact(y)
}

/// Contact `g` at every grid index, evaluated with the current dual `lambda`.
///
/// `refine_socp` deliberately omits the per-time slack, so Algorithm 2 recomputes
/// `g` here for the drop / add / converge logic.
fn contact_on_grid<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    lambda: &Dual,
) -> Vec<f64> {
    (0..grid.len())
        .map(|k| contact_at(cost, grid, gammas, k, lambda))
        .collect()
}

/// Safety backstop for Algorithm 2 (the paper guarantees convergence; the MC
/// targets converge in ≤ 8 iterations). Not a Table III parameter.
const MAX_REFINE_ITERS: usize = 50;

/// Precompute `Γ(t)` over the grid once (`J2Roe` caches nothing — see Design
/// Decision 2). Indexed by grid index.
fn cache_gamma<D: Dynamics>(dynamics: &D, grid: &TimeGrid) -> Vec<SMatrix<f64, N, M>> {
    grid.times().map(|t| dynamics.gamma(t)).collect()
}

/// Solve the fuel-optimal impulsive control problem.
///
/// Wires Algorithm 1 (init) → Algorithm 2 (refine) → Algorithm 3 (extract) with
/// `Γ(t)` caching.
pub fn solve<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
) -> Result<Solution, PlannerError> {
    // --- Input validation (Design Decision 9). ---
    if !grid.dt.is_finite()
        || grid.dt <= 0.0
        || !grid.t_f.is_finite()
        || !grid.t_i.is_finite()
        || grid.t_f <= grid.t_i
    {
        return Err(PlannerError::InvalidInput(
            "grid must satisfy dt > 0 and t_f > t_i (finite)".into(),
        ));
    }
    if params.n_init == 0 || params.n_coarse == 0 {
        return Err(PlannerError::InvalidInput(
            "n_init and n_coarse must be >= 1".into(),
        ));
    }
    let w_norm = w.norm();
    if !w_norm.is_finite() || w_norm <= 0.0 {
        return Err(PlannerError::InvalidInput(
            "target pseudostate w must be nonzero and finite".into(),
        ));
    }

    // --- Cache Γ(t) over the grid. ---
    let gammas = cache_gamma(dynamics, &grid);

    // --- Algorithm 1: initialization (λ_est ∥ w). ---
    let t_est = init::initialize(cost, &grid, &gammas, &w, params);

    // --- Algorithm 2: iterative refinement. ---
    let refined = refine::refine(cost, &grid, &gammas, &w, params, t_est, MAX_REFINE_ITERS)?;

    // --- Algorithm 3: control-input extraction. ---
    let extracted = extract::extract(cost, &grid, &gammas, &w, refined.objective, &refined.t_opt)?;

    Ok(Solution {
        maneuvers: extracted.maneuvers,
        total_dv: extracted.total_dv,
        iterations: refined.iterations,
        residual: extracted.residual,
        lambda: refined.lambda,
    })
}
