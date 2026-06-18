//! Public-API integration tests for the Phase 4 three-step planner.

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::Dynamics;
use koenig_planner::{solve, PlannerError, SolveParams, TimeGrid};
use nalgebra::{SMatrix, SVector};

const N: usize = 6;
const M: usize = 3;

/// Gently rotating, well-conditioned Γ(t) (mirrors the unit-test mock).
struct SpinDyn {
    rate: f64,
}
impl Dynamics for SpinDyn {
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
        let a = self.rate * t;
        let (c, s) = (a.cos(), a.sin());
        SMatrix::<f64, N, M>::from_row_slice(&[
            c,
            -s,
            0.0, //
            s,
            c,
            0.0, //
            0.0,
            0.0,
            1.0, //
            0.5 * c,
            0.0,
            0.5 * s, //
            0.0,
            0.5 * c,
            -0.5 * s, //
            0.5 * s,
            -0.5 * c,
            0.0,
        ])
    }
}

#[test]
fn solve_converges_on_reachable_synthetic_problem() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0); // 61 points
    let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
    let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
    let w = dynamics.gamma(12.0) * ua + dynamics.gamma(47.0) * ub; // reachable
    let cost = Piecewise::new(1.0e12); // Norm2 everywhere
    let params = SolveParams::default();

    let sol = solve(&dynamics, &cost, w, grid, &params).unwrap();

    assert!(sol.residual < 1e-2, "residual {}", sol.residual);
    assert!(sol.iterations >= 1 && sol.iterations <= 50);
    assert!(!sol.maneuvers.is_empty());
    assert!(sol.total_dv > 0.0);
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
}

#[test]
fn solve_rejects_zero_target() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0);
    let w = SVector::<f64, N>::zeros();
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_rejects_degenerate_grid() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 0.0, 0.0); // dt = 0, t_f == t_i
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}
