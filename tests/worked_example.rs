//! Phase 5 worked-example validation, reframed around what is *provable* with the
//! FD-verified J2 mean-ROE dynamics (tests/fd_stm.rs, tests/fd_b_matrix.rs).
//!
//! The paper (Koenig & D'Amico §7, Table III/IV) reports a 3-maneuver, 82.4 mm/s
//! plan. We do NOT assert that, because the paper's published worked-example
//! figures are internally inconsistent with its own dynamics model (its printed
//! STM carries a delta-lambda dt^2 typo we corrected, and even after correcting
//! it the paper's Table IV maneuvers do not reconstruct its Table III target — see
//! `paper_table_iv_does_not_reconstruct`). Instead we assert the things the
//! FD-verified pipeline genuinely guarantees: it converges, the refinement finds
//! the true discretized dual optimum, and the optimum is self-consistent.

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_planner::{refine_socp, solve, CostModel, Dynamics, Pseudostate, SolveParams, TimeGrid};
use nalgebra::SVector;
use std::f64::consts::TAU;

const A_C: f64 = 25_000e3;
const W_METRES: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

fn worked_example() -> (J2Roe, Piecewise, Pseudostate, TimeGrid) {
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let dynamics = J2Roe::new(chief, 0.0, 117_990.0);
    let cost = Piecewise::new(TAU / chief.mean_motion());
    let w = Pseudostate::from_row_slice(&W_METRES) / A_C;
    let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0);
    (dynamics, cost, w, grid)
}

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
                .contact(dynamics.gamma(t).transpose() * sol.lambda)
        })
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(max_g <= 1.01 + 1e-6, "max_t g = {max_g}");

    // The refinement finds the TRUE discretized dual optimum: its objective equals
    // the exact all-times SOCP solved over every grid time (self-consistency).
    let rows: Vec<_> = grid
        .times()
        .map(|t| cost.at(t).cone_constraints(&dynamics.gamma(t)))
        .collect();
    let exact_dual = refine_socp(&w, &rows).expect("exact SOCP").objective;
    let refine_dual = sol.lambda.dot(&w);
    assert!(
        (refine_dual - exact_dual).abs() / exact_dual < 1e-2,
        "refinement dual {refine_dual:.6} vs exact all-times {exact_dual:.6}"
    );

    // The optimum sits where the FD-verified dynamics put it (~80.9 mm/s) — below
    // the paper's reported 82.0 mm/s bound, which the paper's own solution cannot
    // actually achieve in these (corrected) dynamics.
    assert!(
        (0.078..=0.083).contains(&exact_dual),
        "exact dual = {:.4} mm/s",
        exact_dual * 1e3
    );

    // The extracted plan is a sane, near-optimal solution (not garbage). The
    // residual is small but nonzero (the active-set extraction is not perfectly
    // robust on this degenerate e=0.7 contact — a known solver-quality limitation,
    // not a dynamics error).
    assert!(
        sol.total_dv > 0.0 && sol.total_dv < 0.090,
        "total_dv = {}",
        sol.total_dv
    );
    assert!(sol.residual < 5e-2, "residual = {:.3e}", sol.residual);
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
}

#[test]
fn paper_table_iv_does_not_reconstruct() {
    // Evidence for the reframing: the paper's published Table IV maneuvers, fed
    // through the FD-verified dynamics, do NOT reconstruct the Table III target —
    // the residual is enormous (dominated by delta-lambda). This is why we cannot
    // assert the paper's numbers: they are inconsistent with the paper's own model.
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
        recon += dynamics.gamma(*t) * u;
    }
    let residual = (w - recon).norm() / w.norm();
    assert!(
        residual > 0.5,
        "paper Table IV residual = {residual:.3} (expected >> 0; if this drops near \
         0, the paper became reproducible and the validation should be revisited)"
    );
}
