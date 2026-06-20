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

/// Safety backstop for Algorithm 2 (the paper guarantees convergence; in practice
/// the algorithm converges in well under this many iterations). Not a Table III
/// parameter.
const MAX_REFINE_ITERS: usize = 50;

/// Precompute `Γ(t)` over the grid once (`J2Roe` caches nothing internally, so
/// `gamma` is evaluated here and reused). Indexed by grid index.
///
/// # Errors
/// Propagates the first `Dynamics::gamma` failure (out-of-domain chief).
fn cache_gamma<D: Dynamics>(
    dynamics: &D,
    grid: &TimeGrid,
) -> Result<Vec<SMatrix<f64, N, M>>, PlannerError> {
    grid.times().map(|t| dynamics.gamma(t)).collect()
}

/// Shared `solve` preconditions: finite grid with `dt > 0`
/// and `t_f > t_i`, `n_init`/`n_coarse ≥ 1`, and a nonzero finite target `w`.
fn validate_inputs(
    w: &Pseudostate,
    grid: &TimeGrid,
    params: &SolveParams,
) -> Result<(), PlannerError> {
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
    Ok(())
}

/// Run Algorithm 2 (refine) → Algorithm 3 (extract) from a prepared candidate set.
fn run_pipeline<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    params: &SolveParams,
    t_est: Vec<usize>,
) -> Result<Solution, PlannerError> {
    let refined = refine::refine(cost, grid, gammas, w, params, t_est, MAX_REFINE_ITERS)?;
    let extracted = extract::extract(cost, grid, gammas, w, refined.objective, &refined.t_opt)?;
    Ok(Solution {
        maneuvers: extracted.maneuvers,
        total_dv: extracted.total_dv,
        iterations: refined.iterations,
        residual: extracted.residual,
        lambda: refined.lambda,
    })
}

/// Snap arbitrary times to sorted, deduplicated grid indices (nearest grid point,
/// clamped into `[0, len-1]`). Non-finite times are dropped.
fn nearest_grid_indices(grid: &TimeGrid, times: &[f64]) -> Vec<usize> {
    let last = grid.len().saturating_sub(1);
    let mut idx: Vec<usize> = times
        .iter()
        .filter(|t| t.is_finite())
        .map(|&t| ((t - grid.t_i) / grid.dt).round().clamp(0.0, last as f64) as usize)
        .collect();
    idx.sort_unstable();
    idx.dedup();
    idx
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
    validate_inputs(&w, &grid, params)?;
    let gammas = cache_gamma(dynamics, &grid)?;
    // Algorithm 1: initialization (λ_est ∥ w) → the n_init largest-contact times.
    let t_est = init::initialize(cost, &grid, &gammas, &w, params);
    run_pipeline(cost, &grid, &gammas, &w, params, t_est)
}

/// Solve from a caller-supplied initial candidate-time set, **bypassing Algorithm 1**.
///
/// This reproduces the paper's initialization-sensitivity study (Koenig & D'Amico
/// 2020, Fig. 8), whose worst-case `{t_i, t_f}` seed and evenly-spaced seed are
/// deliberately *not* Algorithm-1 (largest-contact) outputs. `initial_times` are
/// snapped to the nearest grid points and deduplicated; refinement (Algorithm 2)
/// then drops/adds times exactly as in [`solve`]. `params.n_init` / `params.n_coarse`
/// are unused on this path (the seed is explicit) but must still be ≥ 1.
///
/// Returns [`PlannerError::InvalidInput`] if no finite initial time lands in range.
pub fn solve_from_initial_times<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
    initial_times: &[f64],
) -> Result<Solution, PlannerError> {
    validate_inputs(&w, &grid, params)?;
    let gammas = cache_gamma(dynamics, &grid)?;
    let t_est = nearest_grid_indices(&grid, initial_times);
    if t_est.is_empty() {
        return Err(PlannerError::InvalidInput(
            "solve_from_initial_times: no finite initial times in range".into(),
        ));
    }
    run_pipeline(cost, &grid, &gammas, &w, params, t_est)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearest_grid_indices_snaps_dedups_and_clamps() {
        let grid = TimeGrid::uniform(0.0, 100.0, 10.0).unwrap(); // indices 0..=10
                                                                 // Endpoints map to the first and last index.
        assert_eq!(nearest_grid_indices(&grid, &[0.0, 100.0]), vec![0, 10]);
        // Nearest rounding, then sorted + deduped: 4->0, 24->2, 26->3.
        assert_eq!(
            nearest_grid_indices(&grid, &[26.0, 4.0, 24.0]),
            vec![0, 2, 3]
        );
        // Out-of-range clamps to the endpoints; non-finite dropped.
        assert_eq!(
            nearest_grid_indices(&grid, &[-50.0, 999.0, f64::NAN, f64::INFINITY]),
            vec![0, 10]
        );
        // Times snapping to the same index collapse to one.
        assert_eq!(nearest_grid_indices(&grid, &[50.0, 51.0, 49.0]), vec![5]);
    }

    #[test]
    fn empty_initial_times_is_rejected() {
        assert!(
            nearest_grid_indices(&TimeGrid::uniform(0.0, 10.0, 1.0).unwrap(), &[f64::NAN])
                .is_empty()
        );
    }
}
