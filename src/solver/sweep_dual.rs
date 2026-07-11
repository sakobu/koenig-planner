//! Batch dual evaluation over a fixed window.
//!
//! [`sweep_dual`] assembles the Γ-cache and conic rows once from
//! `(dynamics, cost, grid)` — all independent of the target `w` — then solves
//! the min-fuel dual ([`refine_socp()`]) for each target. It returns the
//! support-function gauge `c*` (the reachable-set boundary value; \[KD20\] eq. 40)
//! and the dual normal `λ` per target, without rebuilding dynamics or rows per
//! call. Used to trace reachable-set boundaries and Δv cost fields where many
//! targets share one window.

use crate::cost::CostModel;
use crate::dynamics::{Dynamics, J2Roe};
use crate::solver::refine_socp;
use crate::types::{ConicRows, Dual, PlannerError, Pseudostate, TimeGrid};

/// One target's dual solution in a [`sweep_dual`] batch.
#[derive(Debug, Clone)]
pub struct SweepResult {
    /// Support-function gauge `c* = λ*ᵀw` (m/s scale). Meaningful only when
    /// `feasible`; `f64::NAN` for an unreachable target.
    pub c_star: f64,
    /// Optimal dual `λ*` — the outward reachable-set normal. `zeros()` when
    /// infeasible.
    pub lambda: Dual,
    /// `false` when [`refine_socp()`] returned an error — normally an unbounded
    /// dual (the target is unreachable in this window), but also a solver
    /// setup/convergence failure.
    pub feasible: bool,
}

/// Evaluate the min-fuel dual for many targets over one fixed window.
///
/// Assembles the per-time conic rows once (they depend only on `cost` and the
/// Γ(t) from `dynamics`/`grid`, never on the target), then calls [`refine_socp()`]
/// per target. `targets` are dimensionless pseudostates (`δα / a_c`), matching
/// [`crate::solve`]'s `w`.
///
/// # Errors
/// Returns [`PlannerError`] if a Γ(t) evaluation fails for any grid time (a
/// batch-level fault). Per-target unreachability is `SweepResult { feasible:
/// false, .. }`, not an error.
pub fn sweep_dual<C: CostModel>(
    dynamics: &J2Roe,
    cost: &C,
    grid: &TimeGrid,
    targets: &[Pseudostate],
) -> Result<Vec<SweepResult>, PlannerError> {
    // Assemble the conic rows once: rows[k] encodes the cost's sublevel-set
    // constraints on λ from Γ(t_k). Reused for every target ([KD20] eq. 40).
    let mut rows: Vec<ConicRows> = Vec::with_capacity(grid.len());
    for t in grid.times() {
        let gamma_t = dynamics.gamma(t)?;
        rows.push(cost.at(t).cone_constraints(&gamma_t));
    }

    Ok(targets
        .iter()
        .map(|w| match refine_socp(w, &rows) {
            Ok(sol) => SweepResult {
                c_star: sol.objective,
                lambda: sol.lambda,
                feasible: true,
            },
            Err(_) => SweepResult {
                c_star: f64::NAN,
                lambda: Dual::zeros(),
                feasible: false,
            },
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::Piecewise;
    use crate::dynamics::AbsoluteOrbit;
    use approx::assert_relative_eq;
    use std::f64::consts::TAU;

    const A_C: f64 = 25_000e3;

    /// Real J2Roe fixture over one chief period (~64 grid points → fast).
    fn fixture() -> (J2Roe, TimeGrid, Piecewise) {
        let chief = AbsoluteOrbit::new(
            A_C,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        let period = TAU / chief.mean_motion();
        let dynamics = J2Roe::new(chief, 0.0, period).unwrap();
        let grid = TimeGrid::uniform(0.0, period, period / 64.0).unwrap();
        let cost = Piecewise::new(period).unwrap();
        (dynamics, grid, cost)
    }

    // sweep_dual over N targets is byte-identical to N individual refine_socp
    // calls on the same rows — it only hoists the shared assembly.
    #[test]
    fn sweep_matches_individual_refine_socp() {
        let (dynamics, grid, cost) = fixture();
        let w0 = Pseudostate::from_row_slice(&[10.0, 100.0, 10.0, 10.0, 0.0, 10.0]) / A_C;
        let w1 = Pseudostate::from_row_slice(&[5.0, 50.0, 20.0, 0.0, 0.0, 30.0]) / A_C;

        let rows: Vec<ConicRows> = grid
            .times()
            .map(|t| cost.at(t).cone_constraints(&dynamics.gamma(t).unwrap()))
            .collect();
        let r0 = refine_socp(&w0, &rows).unwrap();
        let r1 = refine_socp(&w1, &rows).unwrap();

        let swept = sweep_dual(&dynamics, &cost, &grid, &[w0, w1]).unwrap();

        assert_eq!(swept.len(), 2);
        assert!(swept[0].feasible && swept[1].feasible);
        assert_eq!(swept[0].c_star, r0.objective);
        assert_eq!(swept[1].c_star, r1.objective);
        assert_eq!(swept[0].lambda, r0.lambda);
        assert_eq!(swept[1].lambda, r1.lambda);
    }

    // c* is the support-function gauge → positively homogeneous of degree 1.
    #[test]
    fn c_star_is_positively_homogeneous() {
        let (dynamics, grid, cost) = fixture();
        let w = Pseudostate::from_row_slice(&[10.0, 100.0, 10.0, 10.0, 0.0, 10.0]) / A_C;

        let swept = sweep_dual(&dynamics, &cost, &grid, &[w, w * 2.5]).unwrap();
        assert!(swept[0].feasible && swept[1].feasible);
        assert_relative_eq!(swept[1].c_star, 2.5 * swept[0].c_star, max_relative = 1e-6);
    }

    // A single-candidate-time grid reaches only a 3-D subspace, so a generic
    // 6-D target is unreachable: the dual is unbounded → the point is flagged
    // infeasible rather than erroring the whole batch.
    #[test]
    fn unreachable_target_is_infeasible_not_error() {
        let chief = AbsoluteOrbit::new(
            A_C,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        let dynamics = J2Roe::new(chief, 0.0, 10.0).unwrap();
        // Window (10 s) shorter than dt (100 s) ⇒ a single grid time (t_i only).
        let grid = TimeGrid::uniform(0.0, 10.0, 100.0).unwrap();
        assert_eq!(grid.len(), 1);
        let cost = Piecewise::new(TAU / chief.mean_motion()).unwrap();
        let w = Pseudostate::from_row_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]) / A_C;

        let swept = sweep_dual(&dynamics, &cost, &grid, &[w]).unwrap();
        assert_eq!(swept.len(), 1);
        assert!(!swept[0].feasible);
        assert!(swept[0].c_star.is_nan());
    }
}
