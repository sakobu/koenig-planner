//! Orchestration of the three-step algorithm (Init -> Refine -> Extract).

mod extract;
mod init;
mod refine;

use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{PlannerError, Pseudostate, Solution, SolveParams, TimeGrid};

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
