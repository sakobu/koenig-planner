//! Orchestration of the three-step algorithm (Init -> Refine -> Extract).

mod extract;
mod init;
mod refine;

use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{
    Dual, InvalidInputKind, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N,
};
use nalgebra::{SMatrix, SVector};

/// Contact value `g_{U(1,t)}(Γᵀ(t)·lambda)` at grid index `k`.
///
/// `gammas[k] = Γ(grid.time(k))` is the 6×3 dynamics matrix; the cost methods
/// operate on the control-space projection `Γᵀ(t)·lambda ∈ ℝ³`.
///
/// Ref: \[KD20\] eq. 30 / eq. 27 (contact `g_{U(1,t)}(Γᵀ(t) lambda)`).
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
///
/// Ref: \[KD20\] eq. 30 (max over the candidate grid).
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
        return Err(PlannerError::InvalidInput(InvalidInputKind::Grid {
            t_i: grid.t_i,
            t_f: grid.t_f,
            dt: grid.dt,
        }));
    }
    if params.n_init == 0 || params.n_coarse == 0 {
        return Err(PlannerError::InvalidInput(InvalidInputKind::SolverParams {
            n_init: params.n_init,
            n_coarse: params.n_coarse,
        }));
    }
    let w_norm = w.norm();
    if !w_norm.is_finite() || w_norm <= 0.0 {
        return Err(PlannerError::InvalidInput(InvalidInputKind::Target));
    }
    Ok(())
}

/// Run Algorithm 2 (refine) → Algorithm 3 (extract) from a prepared candidate set.
///
/// Ref: \[KD20\] Algorithm 2 -> Algorithm 3.
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
///
/// Ref: \[KD20\] Algorithms 1/2/3 (full pipeline); eq. 4 (master problem).
///
/// # Errors
/// - [`PlannerError::InvalidInput`] if the grid is non-finite or has `dt <= 0`
///   or `t_f <= t_i`, if `params.n_init` or `params.n_coarse` is `0`, if the
///   target `w` is zero or non-finite, or if the candidate-time set degenerates
///   to empty during refinement.
/// - [`PlannerError::KeplerDivergence`] if evaluating `Γ(t)` runs a Kepler solve
///   whose Newton iteration fails to converge. The built-in `J2Roe` validates
///   its chief elements (all six finite, `a > 0`, `e ∈ [0,1)`), so a malformed
///   chief is rejected as [`PlannerError::InvalidInput`] up front rather than
///   diverging here; a custom `Dynamics` may likewise surface
///   [`PlannerError::InvalidInput`] for an out-of-domain chief.
/// - [`PlannerError::SolverFailed`] if the refinement or min-fuel convex solver
///   fails to set up or solve.
/// - [`PlannerError::NotConverged`] if Algorithm 2 reaches its iteration cap
///   before `max_t g <= 1 + eps_cost`.
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
/// Ref: \[KD20\] Fig. 8 (initialization-sensitivity study; bypasses Algorithm 1).
///
/// # Errors
/// - [`PlannerError::InvalidInput`] if `initial_times` holds no finite value that
///   snaps into the grid.
/// - Every error [`solve`] can return, under the same conditions:
///   [`PlannerError::InvalidInput`] for the shared input validation or an empty
///   refinement set, plus [`PlannerError::KeplerDivergence`],
///   [`PlannerError::SolverFailed`], and [`PlannerError::NotConverged`].
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
            InvalidInputKind::NoInitialTimesInRange,
        ));
    }
    run_pipeline(cost, &grid, &gammas, &w, params, t_est)
}

/// Primer-vector history sampled over a time grid.
///
/// For each grid time `t_k`: `vectors[k] = Γᵀ(t_k)·λ ∈ ℝ³` is the primer vector
/// (\[KD20\] eq. 46 — the dual `λ` mapped into RTN control space), and
/// `magnitudes[k] = g_{U(1,t_k)}(vectors[k])` is its dual-gauge magnitude
/// (dimensionless). This is the paper's Fig. 7 contact curve, paired with the
/// underlying vector.
///
/// The magnitude is `≈ 1` at the optimal maneuver times and `≤ 1 + eps_cost`
/// (Algorithm 2's tolerance) at every candidate time the solve converged over.
/// Evaluated on a grid *denser* than the one solved over it may exceed that
/// bound slightly between solved times, since dual feasibility is only enforced
/// at the solved candidates. Wherever the magnitude touches 1 but no burn is
/// placed, the plan has flexibility.
///
/// `vectors` is the primer, **not** the executed thrust direction in general:
/// the optimal impulse fires along the support image `s_{U(1,t)}(Γᵀλ)`
/// (\[KD20\] eq. 41–42), which is parallel to the primer only for the L2
/// (`Norm2`) gauge; under the `FaceMax` / `Piecewise`-perigee polytope gauge it
/// is a fixed tetrahedral thruster axis, generally not parallel to the primer.
///
/// Ref: \[KD20\] Fig. 7 (contact curve); eq. 30 / 27; eq. 41–42 / 46.
#[derive(Debug, Clone)]
pub struct PrimerHistory {
    /// Sample times `[s]`, one per grid point: `times[k] == grid.time(k)`.
    pub times: Vec<f64>,
    /// Primer vector `p(t_k) = Γᵀ(t_k)·λ`, RTN components `(R, T, N)`.
    pub vectors: Vec<SVector<f64, M>>,
    /// Dual-gauge magnitude `g_{U(1,t_k)}(p(t_k))` (dimensionless).
    pub magnitudes: Vec<f64>,
}

/// Reconstruct the primer-vector history from a converged dual `lambda` over a
/// time grid.
///
/// `lambda` is [`Solution::lambda`] from a [`solve`] call (the eq. 40 dual, the
/// reachable-set normal). The history is evaluated at every `grid` time via the
/// same `Γ(t)` cache and gauge-contact computation Algorithm 2 uses internally,
/// so passing a denser `grid` than the one solved over yields a smoother curve.
/// It evaluates `Γ(t)` once per grid point (`O(grid.len())` dynamics
/// evaluations), so a denser grid trades compute for smoothness.
///
/// Ref: \[KD20\] Fig. 7 (contact curve `g(t)` over the candidate grid); eq. 30.
///
/// # Errors
/// Propagates the first [`Dynamics::gamma`] failure — an out-of-domain chief
/// whose Kepler solve diverges ([`PlannerError::KeplerDivergence`]) or is
/// rejected as non-elliptic ([`PlannerError::InvalidInput`]). Unreachable on the
/// built-in `J2Roe`, which validates its chief in `new`.
pub fn primer_history<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    grid: &TimeGrid,
    lambda: &Dual,
) -> Result<PrimerHistory, PlannerError> {
    let gammas = cache_gamma(dynamics, grid)?;
    let times: Vec<f64> = grid.times().collect();
    let mut vectors = Vec::with_capacity(gammas.len());
    let mut magnitudes = Vec::with_capacity(gammas.len());
    for (g, &t) in gammas.iter().zip(&times) {
        let p = g.transpose() * lambda; // p(t) = Γᵀ(t)·λ ∈ ℝ³, same as `contact_at`
        magnitudes.push(cost.at(t).contact(p));
        vectors.push(p);
    }
    Ok(PrimerHistory {
        times,
        vectors,
        magnitudes,
    })
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
