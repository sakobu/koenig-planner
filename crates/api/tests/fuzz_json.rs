//! Adversarial-JSON tests over the central untrusted-string seam `run_json`.
//! Contract: any input yields Ok or a typed ApiError — never a panic, hang,
//! or OOM. `kind` is BadRequest or Solver (Internal is the HTTP panic-catch
//! kind and must not arise here).

use koenig_damico_planner_api::{run_json, ApiError, ApiErrorKind};
use proptest::prelude::*;

fn cost_spec_json() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(r#"{"type":"norm2"}"#.to_string()),
        Just(r#"{"type":"facemax"}"#.to_string()),
        Just(r#"{"type":"piecewise","period":null,"t_perigee0":null}"#.to_string()),
        Just(r#"{"type":"piecewise","period":39283.0,"t_perigee0":19641.0}"#.to_string()),
    ]
}

prop_compose! {
    // Request-shaped JSON spanning valid AND out-of-domain scalars. Spans/dt
    // are kept small so valid grids stay tiny (fast); extreme-dt DoS is in the
    // corpus (Task 7). Angles are degrees (the api wire convention).
    fn request_shaped_json()(
        a in -1.0e8..1.0e8f64,                 // includes a ≤ 0 (invalid)
        e in -0.5..1.5f64,                     // includes e ≥ 1 and e < 0 (invalid)
        i in -10.0..190.0f64,                  // includes sin i ≈ 0 (invalid)
        raan in -360.0..360.0f64,
        argp in -360.0..360.0f64,
        mean_anom in -360.0..360.0f64,
        t_i in 0.0..100.0f64,
        t_f in 0.0..100.0f64,
        dt in prop_oneof![-50.0..0.0f64, 0.5..50.0f64], // invalid or moderate
        w in proptest::array::uniform6(-1.0e3..1.0e3f64),
        cost in cost_spec_json(),
        init_len in 0usize..8,
    ) -> String {
        let w_str = w.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let init_str = (0..init_len)
            .map(|k| (t_i + k as f64).to_string())
            .collect::<Vec<_>>()
            .join(",");
        format!(
            r#"{{"chief":{{"a":{a},"e":{e},"i":{i},"raan":{raan},"argp":{argp},"mean_anom":{mean_anom}}},"t_i":{t_i},"t_f":{t_f},"dt":{dt},"w_metres":[{w_str}],"cost":{cost},"initial_times":[{init_str}]}}"#
        )
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 256,
        rng_seed: proptest::test_runner::RngSeed::Fixed(0x4B44_3230_5345_4544),
        ..ProptestConfig::default()
    })]

    #[test]
    fn run_json_is_total(s in request_shaped_json()) {
        match run_json(&s) {
            Ok(_) => {}
            Err(ApiError { kind, .. }) => {
                prop_assert!(
                    matches!(kind, ApiErrorKind::BadRequest | ApiErrorKind::Solver),
                    "unexpected kind {:?} for input {}", kind, s
                );
            }
        }
    }
}

/// Concrete malformed inputs must each map to a typed BadRequest, never panic.
#[test]
fn malformed_literals_are_bad_request() {
    let valid_chief = r#"{"a":7000000.0,"e":0.1,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0}"#;
    let cases: Vec<String> = vec![
        "".to_string(),
        "{ not json".to_string(),
        "null".to_string(),
        "[]".to_string(),
        "42".to_string(),
        "{}".to_string(),
        // missing the non-chief fields
        format!(r#"{{"chief":{valid_chief}}}"#),
        // w_metres wrong length (3, needs 6)
        format!(
            r#"{{"chief":{valid_chief},"t_i":0,"t_f":1,"dt":1,"w_metres":[1,2,3],"cost":{{"type":"norm2"}}}}"#
        ),
        // NaN token (serde_json rejects non-finite literals)
        format!(
            r#"{{"chief":{{"a":NaN,"e":0.1,"i":40,"raan":0,"argp":0,"mean_anom":0}},"t_i":0,"t_f":1,"dt":1,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}}}}"#
        ),
        // Infinity token
        format!(
            r#"{{"chief":{valid_chief},"t_i":0,"t_f":Infinity,"dt":1,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}}}}"#
        ),
        // unknown cost tag
        format!(
            r#"{{"chief":{valid_chief},"t_i":0,"t_f":1,"dt":1,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"bogus"}}}}"#
        ),
        // chief wrong type
        format!(
            r#"{{"chief":"x","t_i":0,"t_f":1,"dt":1,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}}}}"#
        ),
        // t_i wrong type
        format!(
            r#"{{"chief":{valid_chief},"t_i":"now","t_f":1,"dt":1,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}}}}"#
        ),
        // truncated
        format!(r#"{{"chief":{valid_chief},"t_i":0,"t_f":1,"dt":1,"#),
    ];
    for c in &cases {
        let err = run_json(c).expect_err(&format!("expected error for: {c}"));
        assert_eq!(err.kind, ApiErrorKind::BadRequest, "input: {c}");
    }
}

/// A request that is valid except for an extra unknown top-level field must be
/// rejected as BadRequest (no silent skip via serde's IgnoredAny path). Also
/// covers an unknown field nested in the `chief` object and in `params`.
#[test]
fn unknown_fields_are_rejected_as_bad_request() {
    let valid_chief = r#"{"a":7000000.0,"e":0.1,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0}"#;
    let cases: Vec<String> = vec![
        // unknown top-level field on SolveRequest
        format!(
            r#"{{"chief":{valid_chief},"t_i":0.0,"t_f":6000.0,"dt":60.0,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}},"junk":1}}"#
        ),
        // unknown field inside the chief (OrbitDto)
        format!(
            r#"{{"chief":{{"a":7000000.0,"e":0.1,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0,"junk":1}},"t_i":0.0,"t_f":6000.0,"dt":60.0,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}}}}"#
        ),
        // unknown field inside params (SolveParamsDto)
        format!(
            r#"{{"chief":{valid_chief},"t_i":0.0,"t_f":6000.0,"dt":60.0,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}},"params":{{"n_coarse":20,"junk":1}}}}"#
        ),
    ];
    for c in &cases {
        let err = run_json(c).expect_err(&format!("expected error for: {c}"));
        assert_eq!(err.kind, ApiErrorKind::BadRequest, "input: {c}");
    }
}

/// Well-formed-but-hostile inputs must not panic, abort, or OOM. They may be
/// Ok or a typed error; the assertion is the absence of a crash.
#[test]
fn stressful_inputs_never_crash() {
    let valid_chief = r#"{"a":7000000.0,"e":0.1,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0}"#;

    // Extreme dt over a wide window → ~1e18 grid points → MAX_GRID_POINTS cap
    // → BadRequest, before any allocation ([KD20] §VIII discretization).
    let extreme_dt = format!(
        r#"{{"chief":{valid_chief},"t_i":0.0,"t_f":1e9,"dt":1e-9,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"facemax"}}}}"#
    );

    // Oversized initial_times (100k entries) over a tiny grid (101 points):
    // they snap+dedup onto the grid, so the SOCP stays bounded (no OOM).
    let big_init = (0..100_000)
        .map(|k| format!("{}.0", k % 100))
        .collect::<Vec<_>>()
        .join(",");
    let oversized = format!(
        r#"{{"chief":{valid_chief},"t_i":0.0,"t_f":6000.0,"dt":60.0,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}},"initial_times":[{big_init}]}}"#
    );

    // Extreme exponents in scalar fields.
    let extreme_vals = r#"{"chief":{"a":1e300,"e":0.999,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0},"t_i":0.0,"t_f":1e9,"dt":1e-9,"w_metres":[1e9,1e9,1e9,1e9,1e9,1e9],"cost":{"type":"facemax"}}"#.to_string();

    for c in [extreme_dt, oversized, extreme_vals] {
        if let Err(e) = run_json(&c) {
            assert!(
                matches!(e.kind, ApiErrorKind::BadRequest | ApiErrorKind::Solver),
                "unexpected kind {:?}",
                e.kind
            );
        }
    }
}

/// An over-large request body is rejected before parsing on the uncapped
/// library path (mirrors the HTTP server's body limit for run_json/py/wasm).
#[test]
fn oversized_request_body_is_rejected() {
    use koenig_damico_planner_api::MAX_REQUEST_BYTES;
    let valid_chief = r#"{"a":7000000.0,"e":0.1,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0}"#;
    // Pad initial_times until the body exceeds the cap.
    let big = (0..)
        .map(|k| (k % 100).to_string())
        .scan(String::new(), |acc, s| {
            if !acc.is_empty() {
                acc.push(',');
            }
            acc.push_str(&s);
            Some(acc.clone())
        })
        .find(|acc| acc.len() > MAX_REQUEST_BYTES)
        .unwrap();
    let body = format!(
        r#"{{"chief":{valid_chief},"t_i":0.0,"t_f":6000.0,"dt":60.0,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}},"initial_times":[{big}]}}"#
    );
    let err = run_json(&body).expect_err("oversized body must reject");
    assert_eq!(err.kind, ApiErrorKind::BadRequest);
}

/// Deeply nested JSON in an UNKNOWN field is now rejected up-front by
/// `deny_unknown_fields` (added on the request DTOs), so the nest is never
/// traversed. Previously the field was silently skipped via serde's IgnoredAny
/// path (iterative/heap-bounded, so no stack overflow, but an unbounded skip
/// on the uncapped library path). The depth-300 nest must reject as BadRequest
/// at the unknown key, with no panic/abort/OOM.
#[test]
fn deeply_nested_unknown_field_is_rejected_without_crash() {
    let valid_chief = r#"{"a":7000000.0,"e":0.1,"i":40.0,"raan":0.0,"argp":0.0,"mean_anom":0.0}"#;
    let depth = 300;
    let deep = format!(
        r#"{{"chief":{valid_chief},"t_i":0.0,"t_f":6000.0,"dt":60.0,"w_metres":[1,2,3,4,5,6],"cost":{{"type":"norm2"}},"junk":{}{}}}"#,
        "[".repeat(depth),
        "]".repeat(depth)
    );
    let err = run_json(&deep).expect_err("unknown field must be rejected");
    assert_eq!(err.kind, ApiErrorKind::BadRequest);
}
