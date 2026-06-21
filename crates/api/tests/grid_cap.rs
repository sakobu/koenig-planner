//! Grid-size cap (audit B1): `run` must reject a request whose grid exceeds
//! `MAX_GRID_POINTS` with `kind == "bad_request"`, before any solve allocation.

use koenig_damico_planner_api::{run, CostSpec, OrbitDto, SolveRequest, MAX_GRID_POINTS};

/// Worked-example chief — a known-valid orbit so `J2Roe::new` accepts it and the
/// request reaches (and trips) the grid-size guard rather than failing earlier.
fn chief() -> OrbitDto {
    OrbitDto {
        a: 25_000e3,
        e: 0.7,
        i: 40.0,
        raan: 358.0,
        argp: 0.0,
        mean_anom: 180.0,
    }
}

/// A tiny `dt` over the valid worked-example window blows the point count far past
/// the cap. `run` must reject it as `bad_request` *before* allocating the Γ-cache
/// (which at this size would be ~17 GB and abort the process).
#[test]
fn oversized_grid_via_tiny_dt_is_bad_request() {
    let req = SolveRequest {
        chief: chief(),
        t_i: 0.0,
        t_f: 117_990.0,
        dt: 0.001, // ~1.18e8 points >> MAX_GRID_POINTS
        w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        cost: CostSpec::Norm2,
        params: None,
        initial_times: None,
    };

    let err = run(req).expect_err("an oversized grid must be rejected");
    assert_eq!(
        err.kind, "bad_request",
        "expected kind=bad_request, got kind={} (message: {})",
        err.kind, err.message
    );
    assert!(
        err.message.contains("grid") && err.message.contains("max"),
        "message should name the grid-size cap, got {:?}",
        err.message
    );
}

/// Exactly one point over the cap: with `t_i=0, dt=1, t_f=MAX_GRID_POINTS`, the
/// endpoint-inclusive point count is `MAX_GRID_POINTS + 1`. DRY against the const;
/// proves the boundary is strict (`>`), not `>=`. Returns before any solve.
#[test]
fn one_over_cap_is_rejected() {
    let req = SolveRequest {
        chief: chief(),
        t_i: 0.0,
        t_f: MAX_GRID_POINTS as f64, // len = MAX_GRID_POINTS + 1
        dt: 1.0,
        w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        cost: CostSpec::Norm2,
        params: None,
        initial_times: None,
    };

    let err = run(req).expect_err("MAX_GRID_POINTS + 1 must be rejected");
    assert_eq!(err.kind, "bad_request", "message: {}", err.message);
    assert!(
        err.message.contains("grid") && err.message.contains("max"),
        "message should name the grid-size cap, got {:?}",
        err.message
    );
}
