use koenig_damico_planner_wasm::version;
use wasm_bindgen_test::*;

#[wasm_bindgen_test]
fn version_is_nonempty() {
    assert!(!version().is_empty());
}

use koenig_damico_planner_wasm::{CostSpec, SolveRequest};

#[wasm_bindgen_test]
fn golden_request_deserializes() {
    let json = r#"{
        "chief": {"a": 25000000.0, "e": 0.7, "i": 40.0, "raan": 358.0, "argp": 0.0, "mean_anom": 180.0},
        "t_i": 0.0, "t_f": 117990.0, "dt": 30.0,
        "w_metres": [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        "cost": {"type": "piecewise"}
    }"#;
    let req: SolveRequest = serde_json::from_str(json).unwrap();
    assert_eq!(req.chief.a, 25_000_000.0);
    assert_eq!(req.w_metres.len(), 6);
    assert!(matches!(
        req.cost,
        CostSpec::Piecewise {
            period: None,
            t_perigee0: None
        }
    ));
}

use koenig_damico_planner_wasm::{solve, solve_json, OrbitDto, SolveOutcome, SolveRequest as Req};
use koenig_damico_planner_wasm::{ApiError, ApiErrorKind, ChiefGeometry, SolveResponse};

fn golden_req() -> Req {
    Req {
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
        cost: koenig_damico_planner_wasm::CostSpec::Piecewise {
            period: None,
            t_perigee0: None,
        },
        params: None,
        initial_times: None,
    }
}

#[wasm_bindgen_test]
fn solve_golden_is_ok_within_bands() {
    match solve(golden_req()) {
        SolveOutcome::Ok { value } => {
            assert!((1..=6).contains(&value.maneuvers.len()));
            assert!(
                value.total_dv > 0.078 && value.total_dv < 0.083,
                "total_dv={}",
                value.total_dv
            );
            assert!(value.residual < 1e-3);
            assert!((1..=50).contains(&value.iterations));
            assert!(value.total_dv.is_finite() && value.lambda.iter().all(|x| x.is_finite()));
            assert_eq!(value.geometry.maneuver_nu.len(), value.maneuvers.len());
            assert!(value.geometry.maneuver_nu.iter().all(|x| x.is_finite()));
            assert!(value.geometry.perigee_window.is_some());
        }
        SolveOutcome::Err { error } => {
            panic!("expected Ok, got err: {:?} {}", error.kind, error.message)
        }
    }
}

#[wasm_bindgen_test]
fn solve_non_elliptic_is_bad_request() {
    let mut req = golden_req();
    req.chief.e = 1.0; // parabolic — rejected upstream
    match solve(req) {
        SolveOutcome::Err { error } => assert_eq!(error.kind, ApiErrorKind::BadRequest),
        SolveOutcome::Ok { .. } => panic!("expected Err for non-elliptic chief"),
    }
}

#[wasm_bindgen_test]
fn solve_outcome_status_tags_are_stable() {
    let ok = SolveOutcome::Ok {
        value: SolveResponse {
            maneuvers: vec![],
            total_dv: 0.0,
            iterations: 0,
            residual: 0.0,
            lambda: [0.0; 6],
            geometry: ChiefGeometry {
                a: 0.0,
                e: 0.0,
                maneuver_nu: vec![],
                perigee_window: None,
            },
        },
    };
    let err = SolveOutcome::Err {
        error: ApiError {
            kind: ApiErrorKind::Solver,
            message: "x".into(),
        },
    };
    assert_eq!(serde_json::to_value(&ok).unwrap()["status"], "ok");
    assert_eq!(serde_json::to_value(&err).unwrap()["status"], "err");
    assert_eq!(
        serde_json::to_value(ApiErrorKind::Internal).unwrap(),
        "internal"
    );
}

#[wasm_bindgen_test]
fn solve_json_roundtrips_and_errors() {
    let json = r#"{"chief":{"a":25000000.0,"e":0.7,"i":40.0,"raan":358.0,"argp":0.0,"mean_anom":180.0},
        "t_i":0.0,"t_f":117990.0,"dt":30.0,"w_metres":[50.0,5000.0,100.0,100.0,0.0,400.0],
        "cost":{"type":"piecewise"}}"#;
    let out = solve_json(json).expect("golden json solves");
    assert!(out.contains("\"maneuvers\""));
    assert!(solve_json("{ not json }").is_err());
}
