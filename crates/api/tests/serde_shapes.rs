//! Serde round-trip and shape tests for the API DTOs.

use koenig_damico_planner_api::{
    run, ApiErrorKind, CostSpec, OrbitDto, SolveParamsDto, SolveRequest,
};

fn minimal_chief() -> OrbitDto {
    OrbitDto {
        a: 25_000e3,
        e: 0.7,
        i: 40.0,
        raan: 358.0,
        argp: 0.0,
        mean_anom: 180.0,
    }
}

fn minimal_request() -> SolveRequest {
    SolveRequest {
        chief: minimal_chief(),
        t_i: 0.0,
        t_f: 117_990.0,
        dt: 30.0,
        w_meters: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        cost: CostSpec::Norm2,
        params: None,
        initial_times: None,
    }
}

/// `SolveRequest` round-trips through JSON unchanged.
#[test]
fn solve_request_round_trips() {
    let req = minimal_request();
    let json = serde_json::to_string(&req).expect("serialize");
    let back: SolveRequest = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(req, back);
}

/// `OrbitDto` round-trips through JSON unchanged.
#[test]
fn orbit_dto_round_trips() {
    let orbit = minimal_chief();
    let json = serde_json::to_string(&orbit).expect("serialize");
    let back: OrbitDto = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(orbit, back);
}

/// `CostSpec::Piecewise` (with `None` fields) round-trips.
#[test]
fn cost_spec_piecewise_round_trips() {
    let cost = CostSpec::Piecewise {
        period: None,
        t_perigee0: None,
    };
    let json = serde_json::to_string(&cost).expect("serialize");
    let back: CostSpec = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(cost, back);
}

/// `CostSpec::Piecewise` with explicit values round-trips.
#[test]
fn cost_spec_piecewise_with_values_round_trips() {
    let cost = CostSpec::Piecewise {
        period: Some(39_338.0),
        t_perigee0: Some(100.0),
    };
    let json = serde_json::to_string(&cost).expect("serialize");
    let back: CostSpec = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(cost, back);
}

/// `SolveParamsDto` with all `None` fields round-trips (default).
#[test]
fn solve_params_dto_defaults_round_trip() {
    let params = SolveParamsDto {
        n_coarse: None,
        n_init: None,
        eps_cost: None,
        eps_remove: None,
    };
    let json = serde_json::to_string(&params).expect("serialize");
    let back: SolveParamsDto = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(params, back);
}

/// The JSON shape of `SolveResponse` has `lambda` as a 6-element array,
/// `maneuvers[i].dv` as a 3-element array, and the three primer-vector arrays
/// present, parallel, and with `primer_rtn` entries of length 3. This is also
/// the only test that runs the standalone `Norm2` cost path, so it doubles as
/// the Norm2-gauge primer regression.
#[test]
fn solve_response_json_shapes() {
    let req = minimal_request(); // CostSpec::Norm2
    let resp = run(req).expect("should solve");
    let json = serde_json::to_string(&resp).expect("serialize");
    let v: serde_json::Value = serde_json::from_str(&json).expect("parse");

    let lambda = v["lambda"].as_array().expect("lambda must be an array");
    assert_eq!(lambda.len(), 6, "lambda must have 6 elements");

    let maneuvers = v["maneuvers"].as_array().expect("maneuvers must be array");
    for (i, m) in maneuvers.iter().enumerate() {
        let dv = m["dv"].as_array().expect("dv must be an array");
        assert_eq!(dv.len(), 3, "maneuver[{i}].dv must have 3 elements");
    }

    // Primer-vector history: three parallel, non-empty arrays; `primer_rtn`
    // entries are 3-vectors; the magnitude reaches the `|p| = 1` bound and never
    // exceeds `1 + eps_cost`.
    let times = v["primer_times"].as_array().expect("primer_times array");
    let mags = v["primer_magnitude"]
        .as_array()
        .expect("primer_magnitude array");
    let rtn = v["primer_rtn"].as_array().expect("primer_rtn array");
    assert!(!times.is_empty(), "primer_times must be non-empty");
    assert_eq!(times.len(), mags.len(), "primer arrays must be parallel");
    assert_eq!(times.len(), rtn.len(), "primer arrays must be parallel");
    assert_eq!(
        rtn[0].as_array().expect("primer_rtn[0] array").len(),
        3,
        "primer_rtn entries must be 3-vectors"
    );
    let max_g = mags
        .iter()
        .map(|g| g.as_f64().expect("magnitude is a number"))
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (1.0 - 1e-3..=1.0 + 0.01 + 1e-6).contains(&max_g),
        "primer max = {max_g} should reach ~1 without exceeding 1 + eps_cost"
    );
}

/// An invalid request (t_f < t_i) maps to `kind == "bad_request"`.
#[test]
fn invalid_time_window_maps_to_bad_request() {
    let req = SolveRequest {
        t_f: 0.0,
        t_i: 100.0, // t_f < t_i — invalid
        ..minimal_request()
    };
    let err = run(req).expect_err("should fail");
    assert_eq!(
        err.kind,
        ApiErrorKind::BadRequest,
        "expected bad_request, got {err}"
    );
}

/// A non-elliptic chief (e >= 1) maps to `kind == "bad_request"`.
#[test]
fn non_elliptic_chief_maps_to_bad_request() {
    let req = SolveRequest {
        chief: OrbitDto {
            e: 1.5, // non-elliptic
            ..minimal_chief()
        },
        ..minimal_request()
    };
    let err = run(req).expect_err("should fail");
    assert_eq!(
        err.kind,
        ApiErrorKind::BadRequest,
        "expected bad_request, got {err}"
    );
}

/// A non-positive semimajor axis maps to `kind == "bad_request"` — NOT `"solver"`.
/// The invalid `a` must be caught at the `J2Roe` gateway, before it poisons the
/// dynamics into a non-finite result that the backstop would misreport as a
/// solver failure.
#[test]
fn nonpositive_semimajor_axis_maps_to_bad_request() {
    let req = SolveRequest {
        chief: OrbitDto {
            a: -25_000e3, // non-positive semimajor axis — physically invalid
            ..minimal_chief()
        },
        ..minimal_request()
    };
    let err = run(req).expect_err("should fail");
    assert_eq!(
        err.kind,
        ApiErrorKind::BadRequest,
        "expected bad_request, got {err}"
    );
}

/// The cost-model wire tags are stable.
#[test]
fn cost_spec_wire_tags_are_stable() {
    assert_eq!(
        serde_json::to_value(CostSpec::Norm2).unwrap()["type"],
        "norm2"
    );
    assert_eq!(
        serde_json::to_value(CostSpec::FaceMax).unwrap()["type"],
        "facemax"
    );
    assert_eq!(
        serde_json::to_value(CostSpec::Piecewise {
            period: None,
            t_perigee0: None
        })
        .unwrap()["type"],
        "piecewise"
    );
}

/// The error-kind wire strings are stable and identical across serde and Display.
#[test]
fn api_error_kind_wire_strings_are_stable() {
    use koenig_damico_planner_api::ApiErrorKind::{BadRequest, Internal, Solver};
    for (kind, want) in [
        (BadRequest, "bad_request"),
        (Solver, "solver"),
        (Internal, "internal"),
    ] {
        assert_eq!(kind.as_str(), want);
        assert_eq!(kind.to_string(), want);
        assert_eq!(
            serde_json::to_value(kind).unwrap(),
            serde_json::Value::String(want.to_string())
        );
    }
}
