//! Public-API integration tests for the three-step planner.

use koenig_damico_planner::cost::Piecewise;
use koenig_damico_planner::dynamics::{AbsoluteOrbit, Dynamics, J2Roe};
use koenig_damico_planner::{solve, solve_from_initial_times, PlannerError, SolveParams, TimeGrid};
use nalgebra::{SMatrix, SVector};
use std::f64::consts::TAU;

const N: usize = 6;
const M: usize = 3;

/// Gently rotating, well-conditioned Γ(t) (mirrors the unit-test mock).
struct SpinDyn {
    rate: f64,
}
impl Dynamics for SpinDyn {
    fn gamma(&self, t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
        let a = self.rate * t;
        let (c, s) = (a.cos(), a.sin());
        Ok(SMatrix::<f64, N, M>::from_row_slice(&[
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
        ]))
    }
}

#[test]
fn solve_converges_on_reachable_synthetic_problem() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap(); // 61 points
    let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
    let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
    let w = dynamics.gamma(12.0).unwrap() * ua + dynamics.gamma(47.0).unwrap() * ub; // reachable
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
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let w = SVector::<f64, N>::zeros();
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_rejects_degenerate_grid() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid {
        t_i: 0.0,
        t_f: 0.0,
        dt: 0.0,
    }; // dt = 0, t_f == t_i
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_rejects_nan_dt() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid {
        t_i: 0.0,
        t_f: 60.0,
        dt: f64::NAN,
    };
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_rejects_infinite_dt() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid {
        t_i: 0.0,
        t_f: 60.0,
        dt: f64::INFINITY,
    };
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_rejects_nan_target() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let w = SVector::<f64, N>::from_row_slice(&[f64::NAN, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

// Ref: [KD20] Algorithm 2; Table III.
#[test]
fn refine_on_real_j2roe_runs_multiple_iterations() {
    // The well-conditioned synthetic problem converges too fast to exercise the
    // drop/add loop body. The real worked-example J2Roe Γ is ill-conditioned
    // (δλ row ~1e3 vs others ~1e-4), so refinement takes several iterations on
    // the degenerate e=0.7 contact. This is a public-API guard that the loop
    // body actually runs on real dynamics.
    const A_C: f64 = 25_000e3;
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let dynamics = J2Roe::new(chief, 0.0, 117_990.0).unwrap();
    let cost = Piecewise::new(TAU / chief.mean_motion());
    let w = SVector::<f64, N>::from_row_slice(&[50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]) / A_C;
    let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0).unwrap();

    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).expect("should solve");
    assert!(
        sol.iterations >= 2,
        "real J2Roe refinement should take >= 2 iterations, got {}",
        sol.iterations
    );
}

// Ref: [KD20] Fig. 8; Table III.
#[test]
fn solve_from_initial_times_endpoints_seed_reconstructs_w() {
    // The paper's worst-case Fig.8 seed (Koenig & D'Amico 2020, p.11) is ONLY the
    // window endpoints {t_i, t_f} — not an Algorithm-1 largest-contact set. It must
    // still converge and reconstruct w via refinement's drop/add.
    const A_C: f64 = 25_000e3;
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let dynamics = J2Roe::new(chief, 0.0, 117_990.0).unwrap();
    let cost = Piecewise::new(TAU / chief.mean_motion());
    let w = SVector::<f64, N>::from_row_slice(&[50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]) / A_C;
    let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0).unwrap();

    let sol = solve_from_initial_times(
        &dynamics,
        &cost,
        w,
        grid,
        &SolveParams::default(),
        &[0.0, 117_990.0],
    )
    .expect("endpoints seed should solve");
    assert!(sol.residual < 1e-3, "residual {}", sol.residual);
    assert!(sol.iterations >= 1 && sol.iterations <= 50);
    assert!(!sol.maneuvers.is_empty());
}

#[test]
fn solve_from_initial_times_rejects_empty_seed() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve_from_initial_times(&dynamics, &cost, w, grid, &SolveParams::default(), &[])
        .unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_propagates_dynamics_gamma_error() {
    // A Dynamics whose gamma always fails must surface through cache_gamma -> solve,
    // proving the fallible-gamma path is wired end-to-end (not swallowed).
    struct FailDyn;
    impl Dynamics for FailDyn {
        fn gamma(&self, _t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
            Err(PlannerError::InvalidInput("boom".into()))
        }
    }
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&FailDyn, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(m) if m == "boom"));
}

#[test]
fn solution_residual_is_unpruned_and_kept_set_residual_is_bounded() {
    // Solution.residual is the residual of the FULL (pre-prune) min-fuel
    // solution. Recomputing it from the *returned* (pruned) maneuvers can yield
    // a larger value, but pruning never decreases it and it stays small. This
    // pins the documented contract.
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
    let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
    let w = dynamics.gamma(12.0).unwrap() * ua + dynamics.gamma(47.0).unwrap() * ub;
    let cost = Piecewise::new(1.0e12);

    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap();

    // Recompute the residual from ONLY the returned (kept) maneuvers.
    let mut w_acc = SVector::<f64, N>::zeros();
    for m in &sol.maneuvers {
        w_acc += dynamics.gamma(m.t).unwrap() * m.dv;
    }
    let kept_residual = (w - w_acc).norm() / w.norm();

    // Reported residual is the small pre-prune one; pruning never decreases it.
    assert!(sol.residual < 1e-2, "reported residual {}", sol.residual);
    assert!(
        kept_residual + 1e-9 >= sol.residual,
        "kept-set residual {kept_residual} should be >= reported {}",
        sol.residual
    );
    assert!(kept_residual < 1e-2, "kept-set residual {kept_residual}");
}

// Ref: [KD20] Algorithm 2.
#[test]
fn solve_reports_not_converged_through_public_api() {
    // eps_cost = -0.5 => convergence target max_t g <= 0.5, which a binding
    // reachable problem (optimal max_t g -> 1.0) can never meet. Refinement
    // exhausts MAX_REFINE_ITERS (50) and surfaces NotConverged through solve().
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
    let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
    let w = dynamics.gamma(12.0).unwrap() * ua + dynamics.gamma(47.0).unwrap() * ub;
    let cost = Piecewise::new(1.0e12);
    let params = SolveParams {
        eps_cost: -0.5,
        ..SolveParams::default()
    };

    let err = solve(&dynamics, &cost, w, grid, &params).unwrap_err();
    match err {
        PlannerError::NotConverged {
            max_iters,
            achieved,
            target,
        } => {
            assert_eq!(max_iters, 50);
            assert!((target - 0.5).abs() < 1e-12, "target {target}");
            assert!(achieved > target, "achieved {achieved} > target {target}");
        }
        other => panic!("expected NotConverged, got {other:?}"),
    }
}

// Ref: [KD20] eq. 40.
#[test]
fn solve_reports_solver_failed_on_unreachable_target() {
    // A constant rank-deficient Gamma (top 3x3 identity; ROE rows 4-6 always
    // zero) leaves lambda_4..6 unconstrained. A target with mass in ROE index 3
    // makes the eq.40 dual (maximize lambda^T w) unbounded => clarabel returns a
    // non-Solved status => SolverFailed propagates through solve().
    struct RankDeficientDyn;
    impl Dynamics for RankDeficientDyn {
        fn gamma(&self, _t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
            let mut g = SMatrix::<f64, N, M>::zeros();
            g[(0, 0)] = 1.0;
            g[(1, 1)] = 1.0;
            g[(2, 2)] = 1.0;
            Ok(g)
        }
    }

    let grid = TimeGrid::uniform(0.0, 60.0, 1.0).unwrap();
    let w = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
    let cost = Piecewise::new(1.0e12);

    let err = solve(&RankDeficientDyn, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(
        matches!(err, PlannerError::SolverFailed(_)),
        "expected SolverFailed, got {err:?}"
    );
}
