//! Batch primal evaluation over a fixed window — the interactive reachable-set
//! engine. Loops the full active-set [`solve()`] (Algorithms 1+2+3) per target
//! and returns cost + dual + confidence + burn count, without the
//! maneuver/geometry payload of a full [`crate::types::Solution`]. Companion to
//! [`crate::solver::sweep_dual()`], which returns the *exact* dual via
//! [`crate::solver::refine_socp()`]; `sweep_solve` returns the (conservative,
//! `<= eps_cost`) primal optimum (\[KD20\] Algorithm 2).

use crate::algorithm::solve;
use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{Dual, Pseudostate, SolveParams, TimeGrid};

/// One target's primal solution in a [`sweep_solve`] batch.
//
// No `#[non_exhaustive]`: the api layer destructures this cross-crate with no
// `..` (the no-drift guard), which `#[non_exhaustive]` would forbid (E0638).
// Matches the existing `SweepResult`.
#[derive(Debug, Clone)]
pub struct SweepSolveResult {
    /// Min-fuel cost `c* = total_dv` (m/s). `f64::NAN` when `!feasible`.
    pub c_star: f64,
    /// Optimal dual `lambda` — the outward reachable-set normal / `grad c*`.
    pub lambda: Dual,
    /// `true` when [`solve()`] returned `Ok` (the target is reachable in this
    /// window). A structurally-unreachable target makes the eq.-40 dual that
    /// [`solve()`] runs internally unbounded, so [`solve()`] returns `Err` —
    /// identical to [`crate::solver::sweep_dual()`]. Reachability is this field,
    /// not a `residual` threshold.
    pub feasible: bool,
    /// Algorithm 2 iteration count (confidence signal).
    pub iterations: u32,
    /// Relative recovery residual `||w_err|| / ||w||` from extract. On a
    /// feasible target this is ~machine-zero (min-fuel enforces `Sum Gamma
    /// Delta-v = w` exactly) — a confidence/accuracy signal, not a
    /// reachability metric. `f64::NAN` when `!feasible`.
    pub residual: f64,
    /// Burn count `maneuvers.len()` (free from extract; drives horizon
    /// burn-count annotations). `0` when `!feasible`.
    pub n_maneuvers: u32,
}

/// Evaluate the min-fuel primal for many targets over one fixed window.
///
/// Loops [`solve()`] per target; a per-target [`crate::types::PlannerError`]
/// becomes `SweepSolveResult { feasible: false, .. }` rather than aborting the
/// batch.
pub fn sweep_solve<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    grid: &TimeGrid,
    targets: &[Pseudostate],
    params: &SolveParams,
) -> Vec<SweepSolveResult> {
    targets
        .iter()
        .map(|w| match solve(dynamics, cost, *w, *grid, params) {
            Ok(s) => SweepSolveResult {
                c_star: s.total_dv,
                lambda: s.lambda,
                feasible: true,
                iterations: s.iterations as u32,
                residual: s.residual,
                n_maneuvers: s.maneuvers.len() as u32,
            },
            Err(_) => SweepSolveResult {
                c_star: f64::NAN,
                lambda: Dual::zeros(),
                feasible: false,
                iterations: 0,
                residual: f64::NAN,
                n_maneuvers: 0,
            },
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::solve;
    use crate::cost::Piecewise;
    use crate::dynamics::{AbsoluteOrbit, J2Roe};
    use std::f64::consts::TAU;

    const A_C: f64 = 25_000e3;

    /// Real J2Roe fixture over one chief period (~64 grid points → fast).
    fn fixture() -> (J2Roe, TimeGrid, Piecewise, SolveParams) {
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
        (dynamics, grid, cost, SolveParams::default())
    }

    #[test]
    fn sweep_solve_matches_individual_solve() {
        let (dynamics, grid, cost, params) = fixture();
        let w0 = Pseudostate::from_row_slice(&[10.0, 100.0, 10.0, 10.0, 0.0, 10.0]) / A_C;
        let w1 = Pseudostate::from_row_slice(&[5.0, 50.0, 20.0, 0.0, 0.0, 30.0]) / A_C;

        let s0 = solve(&dynamics, &cost, w0, grid, &params).unwrap();
        let s1 = solve(&dynamics, &cost, w1, grid, &params).unwrap();

        let swept = sweep_solve(&dynamics, &cost, &grid, &[w0, w1], &params);

        assert_eq!(swept.len(), 2);
        assert!(swept[0].feasible && swept[1].feasible);
        assert_eq!(swept[0].c_star, s0.total_dv);
        assert_eq!(swept[1].c_star, s1.total_dv);
        assert_eq!(swept[0].lambda, s0.lambda);
        assert_eq!(swept[0].n_maneuvers as usize, s0.maneuvers.len());
        assert_eq!(swept[0].iterations as usize, s0.iterations);
        assert_eq!(swept[0].residual, s0.residual);
    }

    #[test]
    fn unreachable_target_is_infeasible() {
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
        let params = SolveParams::default();
        let w = Pseudostate::from_row_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]) / A_C;
        let swept = sweep_solve(&dynamics, &cost, &grid, &[w], &params);
        assert_eq!(swept.len(), 1);
        assert!(!swept[0].feasible);
        assert!(swept[0].c_star.is_nan());
    }
}
