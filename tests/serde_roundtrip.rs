//! Round-trip and shape tests for the optional `serde` feature.
//!
//! Gated on the feature so `cargo test` (no features) ignores this file entirely.
#![cfg(feature = "serde")]

use koenig_damico_planner::{
    dynamics::AbsoluteOrbit, Maneuver, Solution, SolveParams, TimeGrid, M, N,
};
use nalgebra::SVector;

// ---------------------------------------------------------------------------
// AbsoluteOrbit: Serialize + Deserialize — full round-trip
// ---------------------------------------------------------------------------

#[test]
fn absolute_orbit_roundtrips() {
    let orig = AbsoluteOrbit::new(
        25_000e3,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let json = serde_json::to_string(&orig).expect("serialize AbsoluteOrbit");
    let back: AbsoluteOrbit = serde_json::from_str(&json).expect("deserialize AbsoluteOrbit");
    assert_eq!(orig, back);
}

// ---------------------------------------------------------------------------
// TimeGrid: Serialize + Deserialize — full round-trip
// ---------------------------------------------------------------------------

#[test]
fn time_grid_roundtrips() {
    let orig = TimeGrid::uniform(0.0, 117990.0, 30.0).unwrap();
    let json = serde_json::to_string(&orig).expect("serialize TimeGrid");
    let back: TimeGrid = serde_json::from_str(&json).expect("deserialize TimeGrid");
    assert_eq!(orig, back);
}

// ---------------------------------------------------------------------------
// SolveParams: Serialize + Deserialize — full round-trip
// ---------------------------------------------------------------------------

#[test]
fn solve_params_default_roundtrips() {
    let orig = SolveParams::default();
    let json = serde_json::to_string(&orig).expect("serialize SolveParams");
    let back: SolveParams = serde_json::from_str(&json).expect("deserialize SolveParams");
    assert_eq!(orig.n_coarse, back.n_coarse);
    assert_eq!(orig.n_init, back.n_init);
    assert_eq!(orig.eps_cost, back.eps_cost);
    assert_eq!(orig.eps_remove, back.eps_remove);
}

// ---------------------------------------------------------------------------
// Maneuver: Serialize-only — assert `dv` is a flat 3-element JSON array
// ---------------------------------------------------------------------------

#[test]
fn maneuver_dv_serializes_as_flat_array() {
    let m = Maneuver {
        t: 1234.5,
        dv: SVector::<f64, M>::new(0.1, 0.2, 0.3),
    };
    let json = serde_json::to_string(&m).expect("serialize Maneuver");
    let v: serde_json::Value = serde_json::from_str(&json).expect("parse Maneuver JSON");

    // `dv` must be a 3-element array (flat, from nalgebra/serde-serialize)
    let dv = v.get("dv").expect("dv key present");
    assert!(dv.is_array(), "dv must be a JSON array");
    let arr = dv.as_array().unwrap();
    assert_eq!(arr.len(), M, "dv must have {M} elements");
    assert_eq!(arr[0].as_f64().unwrap(), 0.1);
    assert_eq!(arr[1].as_f64().unwrap(), 0.2);
    assert_eq!(arr[2].as_f64().unwrap(), 0.3);
}

// ---------------------------------------------------------------------------
// Solution: Serialize-only — assert `lambda` is a flat 6-element JSON array
// ---------------------------------------------------------------------------

#[test]
fn solution_lambda_serializes_as_flat_array() {
    let lambda = SVector::<f64, N>::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    let sol = Solution {
        maneuvers: vec![],
        total_dv: 0.5,
        iterations: 3,
        residual: 1e-14,
        lambda,
    };
    let json = serde_json::to_string(&sol).expect("serialize Solution");
    let v: serde_json::Value = serde_json::from_str(&json).expect("parse Solution JSON");

    // `lambda` must be a 6-element array (flat, from nalgebra/serde-serialize)
    let lam = v.get("lambda").expect("lambda key present");
    assert!(lam.is_array(), "lambda must be a JSON array");
    let arr = lam.as_array().unwrap();
    assert_eq!(arr.len(), N, "lambda must have {N} elements");
    for (idx, expected) in (1..=6_u64).enumerate() {
        assert_eq!(arr[idx].as_f64().unwrap(), expected as f64);
    }
}
