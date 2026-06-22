//! Worked-example validation, scoped to what is *provable* with the
//! FD-verified J2 mean-ROE dynamics (tests/fd_stm.rs, tests/fd_b_matrix.rs).
//!
//! The paper (Koenig & D'Amico §7, Table III/IV) reports a 3-maneuver, 82.4 mm/s
//! plan. This test does not assert those figures, because the paper's published
//! worked-example figures are internally inconsistent with its own dynamics model
//! (its printed STM carries a delta-lambda dt^2 typo — this crate uses the
//! corrected linear term — and even after correcting it the paper's Table IV
//! maneuvers do not reconstruct its Table III target; see
//! `paper_table_iv_does_not_reconstruct`). This test instead asserts the things
//! the finite-difference-verified pipeline guarantees: it converges, the
//! refinement finds the true discretized dual optimum, and the optimum is
//! self-consistent.

use approx::{assert_abs_diff_eq, assert_relative_eq};
use koenig_damico_planner::cost::Piecewise;
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::{
    refine_socp, solve, CostModel, Dynamics, Pseudostate, SolveParams, TimeGrid,
};
use nalgebra::SVector;
use std::f64::consts::TAU;

const A_C: f64 = 25_000e3;
const W_METRES: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

// Ref: [KD20] Table III; eq. 49.
fn worked_example() -> (J2Roe, Piecewise, Pseudostate, TimeGrid) {
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let dynamics = J2Roe::new(chief, 0.0, 117_990.0).unwrap();
    let cost = Piecewise::new(TAU / chief.mean_motion()).unwrap();
    let w = Pseudostate::from_row_slice(&W_METRES) / A_C;
    let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0).unwrap();
    (dynamics, cost, w, grid)
}

// Ref: [KD20] Table III; eq. 40; eq. 30 / 27; eq. 42.
#[test]
fn worked_example_is_self_consistent() {
    let (dynamics, cost, w, grid) = worked_example();
    assert_eq!(grid.len(), 3934, "Table III grid should have 3934 times");
    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).expect("should solve");

    // Converges within the iteration cap (paper: ~3 iterations).
    assert!(
        (1..=50).contains(&sol.iterations),
        "iterations = {}",
        sol.iterations
    );

    // The optimality criterion is satisfied: max_t g <= 1 + eps_cost.
    let max_g = grid
        .times()
        .map(|t| {
            cost.at(t)
                .contact(dynamics.gamma(t).unwrap().transpose() * sol.lambda)
        })
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(max_g <= 1.01 + 1e-6, "max_t g = {max_g}");

    // The refinement finds the TRUE discretized dual optimum: its objective equals
    // the exact all-times SOCP solved over every grid time (self-consistency).
    let rows: Vec<_> = grid
        .times()
        .map(|t| cost.at(t).cone_constraints(&dynamics.gamma(t).unwrap()))
        .collect();
    let exact_dual = refine_socp(&w, &rows).expect("exact SOCP").objective;
    assert!(
        exact_dual > 0.0,
        "exact_dual must be positive, got {exact_dual}"
    );
    let refine_dual = sol.lambda.dot(&w);
    assert!(
        (refine_dual - exact_dual).abs() / exact_dual < 1e-2,
        "refinement dual {refine_dual:.6} vs exact all-times {exact_dual:.6}"
    );

    // The optimum (~80.9 mm/s) sits below the paper's reported 82.0 mm/s bound under
    // the corrected dynamics used here.
    assert!(
        (0.078..=0.083).contains(&exact_dual),
        "exact dual = {:.4} mm/s",
        exact_dual * 1e3
    );

    // The direct min-fuel SOCP recovers w to ~0 residual with a small maneuver set.
    // `total_dv` is the min-fuel primal objective (the paper's "delta-v cost" c*):
    // ~81.4 mm/s here, sitting just above the 80.9 mm/s dual lower bound (the ~0.6%
    // gap is the eps_cost = 1% convergence tolerance). Both lie in the [78, 83] band.
    assert!(
        sol.total_dv > 0.078 && sol.total_dv < 0.083,
        "total_dv = {} (expected ~81.4 mm/s, the min-fuel primal objective)",
        sol.total_dv
    );
    assert!(
        sol.residual < 1e-3,
        "residual = {:.3e} (target: << 0.1%)",
        sol.residual
    );
    // <= N = 6: the R^6 dual needs at most N active contacts at the optimum
    // (Caratheodory), so >6 post-prune maneuvers signals a pruning or
    // sparsity regression. Observed here: 3.
    assert!(
        sol.maneuvers.len() <= 6,
        "expected a sparse maneuver set (<= 6), got {}",
        sol.maneuvers.len()
    );
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
    // --- Tight regression pin: the converged worked-example solution ([KD20] Table III/IV). ---
    // The coarse bands above are a paper-bound sanity check; these catch silent science drift
    // (a 1–3% Δv shift or a residual blow-up) that the bands would pass.
    assert_relative_eq!(sol.total_dv, 8.135_423_473_009_615e-2, max_relative = 1e-9);
    assert!(
        sol.residual < 1e-10,
        "residual {:.3e} regressed above the 1e-10 ceiling (converges to ~1.1e-14)",
        sol.residual
    );
    assert_eq!(
        sol.maneuvers.len(),
        3,
        "maneuver count regressed from the expected 3"
    );
    let mut man: Vec<&_> = sol.maneuvers.iter().collect();
    man.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    let expected: [(f64, f64); 3] = [
        (14_460.0, 3.952_133_923_056_219e-2),
        (101_970.0, 2.487_437_546_693_282e-2),
        (104_190.0, 1.695_846_939_211_844_6e-2),
    ];
    for (m, (t, dv)) in man.iter().zip(expected) {
        assert_abs_diff_eq!(m.t, t, epsilon = 1e-6);
        assert_relative_eq!(m.dv.norm(), dv, max_relative = 1e-7);
    }
}

// Ref: [H25] Table 2/3 / eq. 69; eq. 19; [KD20] eq. 40.
#[test]
fn hunter_l2_cross_check_recovers_w() {
    // Hunter & D'Amico 2025 "Sequential Formulation Validation": identical J2 ROE
    // dynamics, pure L2 cost. The dual lower bound is correct (~2.48e-4 m/s in the
    // dynamics used here; the paper's 2.294e-4 is not reproduced). This test asserts
    // that the min-fuel extractor recovers w to <0.01% residual at the
    // self-consistent dual.
    let e_x: f64 = -0.658;
    let e_y: f64 = -0.239;
    let e = (e_x * e_x + e_y * e_y).sqrt();
    let argp = e_y.atan2(e_x); // atan2 -> -2.7932 rad (-160 deg = 200 deg mod 360)
    let u0 = 65.0_f64.to_radians();
    let mean_anom = u0 - argp; // u0 = argp + M -> +3.9277 rad (= -135 deg mod 360); Kepler wraps it
    let chief = AbsoluteOrbit::new(
        A_C,
        e,
        51.0_f64.to_radians(),
        30.0_f64.to_radians(),
        argp,
        mean_anom,
    );
    let dynamics = J2Roe::new(chief, 0.0, 39_000.0).unwrap();
    let cost = Piecewise::new(1.0e18).unwrap(); // pure Norm2 (no perigee window ever active)
    let w = Pseudostate::from_row_slice(&[0.66, -1.52, -0.38, -1.44, 0.29, -0.91]) / A_C;
    let grid = TimeGrid::uniform(0.0, 39_000.0, 10.0).unwrap();
    assert_eq!(grid.len(), 3901);

    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).expect("should solve");

    // Self-consistency: refinement objective equals the exact all-times dual.
    let rows: Vec<_> = grid
        .times()
        .map(|t| cost.at(t).cone_constraints(&dynamics.gamma(t).unwrap()))
        .collect();
    let exact_dual = refine_socp(&w, &rows).expect("exact SOCP").objective;
    assert!(
        exact_dual > 0.0,
        "exact_dual must be positive, got {exact_dual}"
    );
    assert!(
        (sol.lambda.dot(&w) - exact_dual).abs() / exact_dual < 1e-2,
        "refine dual {} vs exact {exact_dual}",
        sol.lambda.dot(&w)
    );

    // Extraction reconstructs w to < 0.01% residual.
    assert!(sol.residual < 1e-4, "residual = {:.3e}", sol.residual);
    // Optimum under these dynamics (~2.487e-4 m/s), not the paper's 2.294e-4 bound.
    assert!(
        (2.4e-4..=2.6e-4).contains(&sol.total_dv),
        "total_dv = {:.4e} m/s",
        sol.total_dv
    );
    assert!((1..=50).contains(&sol.iterations));
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
    // --- Tight regression pin: the Hunter L2 cross-check converged solution. ---
    assert_relative_eq!(sol.total_dv, 2.487_623_042_218_083e-4, max_relative = 1e-9);
    assert!(
        sol.residual < 1e-7,
        "residual {:.3e} regressed above the 1e-7 ceiling (converges to ~4.6e-9)",
        sol.residual
    );
    assert_eq!(
        sol.maneuvers.len(),
        5,
        "maneuver count regressed from the expected 5"
    );
    let mut man: Vec<&_> = sol.maneuvers.iter().collect();
    man.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    let expected: [(f64, f64); 5] = [
        (0.0, 8.394_709_997_324_14e-5),
        (6_380.0, 2.487_701_503_472_918e-7),
        (7_170.0, 1.653_549_728_822_607e-5),
        (18_460.0, 1.456_169_138_713_659_2e-4),
        (20_530.0, 2.350_619_843_350_491_7e-6),
    ];
    for (m, (t, dv)) in man.iter().zip(expected) {
        assert_abs_diff_eq!(m.t, t, epsilon = 1e-6);
        assert_relative_eq!(m.dv.norm(), dv, max_relative = 1e-7);
    }
}

// Ref: [KD20] Table IV vs Table III; eq. 11.
#[test]
fn paper_table_iv_does_not_reconstruct() {
    // Demonstrates that the paper's published Table IV maneuvers, fed
    // through the FD-verified dynamics, do NOT reconstruct the Table III target —
    // the residual is enormous (dominated by delta-lambda), so the paper's published
    // figures are not reproduced by the paper's own model.
    let (dynamics, _cost, w, _grid) = worked_example();
    let table_iv = [
        (
            16_050.0,
            SVector::<f64, 3>::new(9.68e-3, -23.02e-3, -25.56e-3),
        ),
        (
            23_280.0,
            SVector::<f64, 3>::new(0.00e-3, -0.40e-3, -0.04e-3),
        ),
        (
            107_100.0,
            SVector::<f64, 3>::new(16.51e-3, 15.68e-3, 40.26e-3),
        ),
    ];
    let mut recon = SVector::<f64, 6>::zeros();
    for (t, u) in &table_iv {
        recon += dynamics.gamma(*t).unwrap() * u;
    }
    let residual = (w - recon).norm() / w.norm();
    assert!(
        residual > 0.5,
        "paper Table IV residual = {residual:.3} (expected >> 0; the paper's figures are not reproduced by its own model)"
    );
}
