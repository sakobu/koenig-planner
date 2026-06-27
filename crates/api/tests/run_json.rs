use koenig_damico_planner_api::{run_json, ApiErrorKind};
use serde_json::Value;

const WORKED_EXAMPLE_JSON: &str = r#"{
    "chief": {"a": 25000000.0, "e": 0.7, "i": 40.0, "raan": 358.0, "argp": 0.0, "mean_anom": 180.0},
    "t_i": 0.0,
    "t_f": 117990.0,
    "dt": 30.0,
    "w_meters": [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
    "cost": {"type": "piecewise"}
}"#;

#[test]
fn run_json_golden_roundtrip() {
    let output = run_json(WORKED_EXAMPLE_JSON).expect("worked example should solve");
    let v: Value = serde_json::from_str(&output).expect("output must be valid JSON");

    let maneuvers = v["maneuvers"]
        .as_array()
        .expect("maneuvers must be an array");
    assert!(
        (1..=6).contains(&maneuvers.len()),
        "expected 1..=6 maneuvers, got {}",
        maneuvers.len()
    );

    let total_dv = v["total_dv"].as_f64().expect("total_dv must be a number");
    assert!(
        (0.078..0.083).contains(&total_dv),
        "total_dv {total_dv} out of range (0.078, 0.083)"
    );

    let residual = v["residual"].as_f64().expect("residual must be a number");
    assert!(residual < 1e-3, "residual {residual} must be < 1e-3");

    let lambda = v["lambda"].as_array().expect("lambda must be an array");
    assert_eq!(lambda.len(), 6, "lambda must have 6 elements");
}

#[test]
fn run_json_malformed_returns_bad_request() {
    let err = run_json("{ not json").unwrap_err();
    assert_eq!(err.kind, ApiErrorKind::BadRequest);
}

// The WASM/HTTP frontends return `ApiError` directly as a JSON error body, so it
// must serialize to the `{kind, message}` wire shape.
#[test]
fn api_error_serializes_to_wire_json() {
    let err = run_json("{ not json").unwrap_err();
    let json = serde_json::to_string(&err).expect("ApiError must serialize");
    let v: Value = serde_json::from_str(&json).expect("serialized ApiError must be valid JSON");
    assert_eq!(v["kind"], "bad_request");
    assert!(
        v["message"]
            .as_str()
            .expect("message must be a string")
            .contains("invalid request JSON"),
        "message should describe the parse failure, got {:?}",
        v["message"]
    );
}
