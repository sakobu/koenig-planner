//! Primer-vector history: the dual-gauge contact curve
//! `g_{U(1,t)}(Γᵀ(t)·λ)` — the paper's Fig. 7 — exposed as a per-solve output,
//! together with the primer vector `p(t) = Γᵀ(t)·λ` itself (RTN).
//!
//! The magnitude is `<= 1 + eps_cost` everywhere (Algorithm 2 convergence) and
//! `~= 1` at the optimal maneuver times (complementary slackness): wherever the
//! curve touches 1 but no burn is placed, the plan has flexibility.
//!
//! Ref: \[KD20\] Fig. 7 (contact curve); eq. 30 / 27; eq. 40.

use approx::assert_relative_eq;
use koenig_damico_planner::cost::{FaceMax, Piecewise};
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::{
    primer_history, solve, CostModel, Dynamics, Pseudostate, SolveParams, SublevelSet, TimeGrid,
};
use std::f64::consts::TAU;

/// Constant `FaceMax` cost model (FaceMax at every time) — mirrors the private
/// adapter the api uses, so a pure-FaceMax solve can be driven from a test.
struct ConstFaceMax(FaceMax);
impl CostModel for ConstFaceMax {
    fn at(&self, _t: f64) -> &dyn SublevelSet {
        &self.0
    }
}

const A_C: f64 = 25_000e3;
const W_METRES: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

// Ref: [KD20] Table III; eq. 49 (identical to tests/worked_example.rs).
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

// The history is aligned to the grid, finite, dual-feasible, and the reported
// magnitude is exactly the gauge contact of the reported (RTN) primer vector.
#[test]
fn primer_history_is_grid_aligned_and_dual_feasible() {
    let (dynamics, cost, w, grid) = worked_example();
    let params = SolveParams::default();
    let sol = solve(&dynamics, &cost, w, grid, &params).expect("worked example should solve");

    let hist = primer_history(&dynamics, &cost, &grid, &sol.lambda).expect("primer history");

    // One sample per grid time, in grid order.
    assert_eq!(hist.times.len(), grid.len());
    assert_eq!(hist.magnitudes.len(), grid.len());
    assert_eq!(hist.vectors.len(), grid.len());
    for (k, t) in grid.times().enumerate() {
        assert_eq!(hist.times[k], t, "times[{k}] should equal grid.time({k})");
    }

    // Everything finite (the api finite-guard relies on this).
    assert!(hist.magnitudes.iter().all(|m| m.is_finite()));
    assert!(hist.vectors.iter().all(|v| v.iter().all(|x| x.is_finite())));

    // The reported vector is the primer `p(t) = Γᵀ(t)·λ`, and the reported
    // magnitude is exactly the active gauge's contact of it.
    for k in 0..grid.len() {
        let t = grid.time(k);
        let p = dynamics.gamma(t).unwrap().transpose() * sol.lambda;
        assert_relative_eq!((hist.vectors[k] - p).norm(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(
            hist.magnitudes[k],
            cost.at(t).contact(p),
            max_relative = 1e-12
        );
    }

    // Dual feasibility: g <= 1 + eps_cost everywhere (Algorithm 2 stop criterion)
    // and the curve actually reaches ~1 (not degenerately slack).
    let max_g = hist
        .magnitudes
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    assert!(
        max_g <= 1.0 + params.eps_cost + 1e-6,
        "max_t g = {max_g} exceeds 1 + eps_cost"
    );
    assert!(max_g >= 1.0 - 1e-3, "max_t g = {max_g} should reach ~1");
}

// Complementary slackness: every optimal maneuver sits at a near-active contact,
// so its primer magnitude is within eps_cost of 1. This pins that we expose the
// quantity whose "touches 1" structure reveals plan flexibility.
#[test]
fn primer_touches_one_at_each_maneuver() {
    let (dynamics, cost, w, grid) = worked_example();
    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).expect("should solve");
    let hist = primer_history(&dynamics, &cost, &grid, &sol.lambda).expect("primer history");

    assert!(!sol.maneuvers.is_empty(), "worked example has maneuvers");
    for m in &sol.maneuvers {
        let k = hist
            .times
            .iter()
            .position(|&t| (t - m.t).abs() < 1e-6)
            .unwrap_or_else(|| panic!("maneuver time {} not on the grid", m.t));
        let g = hist.magnitudes[k];
        assert!(
            (0.98..=1.0 + 2e-2).contains(&g),
            "maneuver at t={} has primer magnitude {g} (expected ~1)",
            m.t
        );
    }
}

// The worked-example maneuvers all land on Norm2 (non-perigee) times, so the
// Piecewise tests above never exercise an *active* FaceMax sample. A pure-FaceMax
// cost makes every maneuver time a FaceMax-active contact, pinning the primer's
// touches-1 / dual-feasibility structure through the polytope gauge too.
#[test]
fn primer_touches_one_at_each_maneuver_under_facemax() {
    let (dynamics, _piecewise, w, grid) = worked_example();
    let cost = ConstFaceMax(FaceMax);
    let params = SolveParams::default();
    let sol = solve(&dynamics, &cost, w, grid, &params).expect("FaceMax should solve");
    let hist = primer_history(&dynamics, &cost, &grid, &sol.lambda).expect("primer history");

    assert!(hist.magnitudes.iter().all(|m| m.is_finite()));
    let max_g = hist
        .magnitudes
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    assert!(
        max_g <= 1.0 + params.eps_cost + 1e-6,
        "FaceMax max_t g = {max_g} exceeds 1 + eps_cost"
    );
    assert!(
        max_g >= 1.0 - 1e-3,
        "FaceMax max_t g = {max_g} should reach ~1"
    );

    assert!(!sol.maneuvers.is_empty(), "FaceMax solve has maneuvers");
    for m in &sol.maneuvers {
        let k = hist
            .times
            .iter()
            .position(|&t| (t - m.t).abs() < 1e-6)
            .unwrap_or_else(|| panic!("maneuver time {} not on the grid", m.t));
        let g = hist.magnitudes[k];
        assert!(
            (0.98..=1.0 + 2e-2).contains(&g),
            "FaceMax maneuver at t={} has primer magnitude {g} (expected ~1)",
            m.t
        );
    }
}
