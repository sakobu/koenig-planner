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

/// Solve the fuel-optimal impulsive control problem.
///
/// Wires Algorithm 1 (init) -> Algorithm 2 (refine) -> Algorithm 3 (extract)
/// with `Gamma(t)` caching. Implemented in Phases 4-5.
#[allow(unused_variables)]
pub fn solve<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
) -> Result<Solution, PlannerError> {
    unimplemented!("Phases 4-5 wire init -> refine -> extract")
}
