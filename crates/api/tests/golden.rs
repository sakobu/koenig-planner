//! Golden worked-example test.
//!
//! Asserts the same bands as `tests/worked_example.rs` in the core crate —
//! these are FD-verified, solver-dependent ranges rather than exact bytes.

use koenig_planner_api::{run, CostSpec, OrbitDto, SolveParamsDto, SolveRequest};

/// Canonical Koenig & D'Amico (2020) worked example (Table III).
///
/// Chief: a=25 000 km, e=0.7, i=40°, Ω=358°, ω=0°, M₀=180° (apogee).
/// Control window: [0, 117 990] s, dt = 30 s.
/// Target pseudostate w [m]: [50, 5000, 100, 100, 0, 400].
/// Cost: Piecewise with period derived from the chief (None → TAU/n).
#[test]
fn golden_worked_example() {
    let req = SolveRequest {
        chief: OrbitDto {
            a: 25_000e3,
            e: 0.7,
            i: 40.0,
            raan: 358.0,
            argp: 0.0,
            mean_anom: 180.0,
        },
        t_i: 0.0,
        t_f: 117_990.0,
        dt: 30.0,
        w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        cost: CostSpec::Piecewise {
            period: None,
            t_perigee0: None,
        },
        params: None,
        initial_times: None,
    };

    let resp = run(req).expect("worked example should solve");

    // Carathéodory bound: at most 6 maneuvers (N = 6 ROE dimensions).
    assert!(
        (1..=6).contains(&resp.maneuvers.len()),
        "expected 1–6 maneuvers, got {}",
        resp.maneuvers.len()
    );

    // FD-verified total Δv band [m/s] (observed ≈ 0.0808).
    assert!(
        resp.total_dv > 0.078 && resp.total_dv < 0.083,
        "total_dv = {} is outside [0.078, 0.083]",
        resp.total_dv
    );

    // Relative residual well below 1 % (observed ≈ 1e-14).
    assert!(
        resp.residual < 1e-3,
        "residual = {} is too large (expected < 1e-3)",
        resp.residual
    );

    // Algorithm 2 should converge within the 50-iteration backstop.
    assert!(
        (1..=50).contains(&resp.iterations),
        "iterations = {} is outside [1, 50]",
        resp.iterations
    );

    // All outputs must be finite (the finite-guard in `run` ensures this, but
    // asserting here keeps the test self-contained).
    assert!(
        resp.lambda.iter().all(|x| x.is_finite()),
        "lambda contains non-finite values: {:?}",
        resp.lambda
    );
    for (i, m) in resp.maneuvers.iter().enumerate() {
        assert!(
            m.dv.iter().all(|x| x.is_finite()),
            "maneuver {} dv contains non-finite values: {:?}",
            i,
            m.dv
        );
    }
}

/// FaceMax cost model dispatch — exercises the `ConstFaceMax` adapter via
/// `run()`.  Reuses the worked-example chief/grid/`w_metres`.  FaceMax is a
/// different cost from Piecewise so no specific Δv band is asserted; only
/// finiteness, sparsity (1–6 maneuvers), and positivity are checked.
#[test]
fn facemax_run_ok() {
    let req = SolveRequest {
        chief: OrbitDto {
            a: 25_000e3,
            e: 0.7,
            i: 40.0,
            raan: 358.0,
            argp: 0.0,
            mean_anom: 180.0,
        },
        t_i: 0.0,
        t_f: 117_990.0,
        dt: 30.0,
        w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        cost: CostSpec::FaceMax,
        params: None,
        initial_times: None,
    };

    let resp = run(req).expect("FaceMax on worked-example orbit should solve");

    // Carathéodory bound: at most 6 maneuvers (N = 6 ROE dimensions).
    assert!(
        (1..=6).contains(&resp.maneuvers.len()),
        "expected 1–6 maneuvers, got {}",
        resp.maneuvers.len()
    );

    // Total Δv must be strictly positive.
    assert!(
        resp.total_dv.is_finite() && resp.total_dv > 0.0,
        "total_dv = {} is not finite and positive",
        resp.total_dv
    );

    // Residual must be finite.
    assert!(
        resp.residual.is_finite(),
        "residual = {} is not finite",
        resp.residual
    );

    // Every dv component and every lambda component must be finite.
    for (i, m) in resp.maneuvers.iter().enumerate() {
        assert!(
            m.dv.iter().all(|x| x.is_finite()),
            "maneuver {} dv contains non-finite values: {:?}",
            i,
            m.dv
        );
    }
    assert!(
        resp.lambda.iter().all(|x| x.is_finite()),
        "lambda contains non-finite values: {:?}",
        resp.lambda
    );
}

/// A `n_coarse: 0` request must map to `kind == "bad_request"` (Fix 1:
/// `PlannerError::InvalidInput` routed correctly).
#[test]
fn n_coarse_zero_is_bad_request() {
    let req = SolveRequest {
        chief: OrbitDto {
            a: 25_000e3,
            e: 0.7,
            i: 40.0,
            raan: 358.0,
            argp: 0.0,
            mean_anom: 180.0,
        },
        t_i: 0.0,
        t_f: 117_990.0,
        dt: 30.0,
        w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        cost: CostSpec::Piecewise {
            period: None,
            t_perigee0: None,
        },
        params: Some(SolveParamsDto {
            n_coarse: Some(0),
            n_init: None,
            eps_cost: None,
            eps_remove: None,
        }),
        initial_times: None,
    };

    let err = run(req).expect_err("n_coarse=0 should return an error");
    assert_eq!(
        err.kind, "bad_request",
        "expected kind=bad_request, got kind={} (message: {})",
        err.kind, err.message
    );
}
