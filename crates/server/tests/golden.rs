//! HTTP-level contract tests, driven in-process via `tower::ServiceExt::oneshot`
//! (no socket). Reuses the Koenig & D'Amico (2020) worked-example golden fixture.

use axum::body::Body;
use axum::http::{header, Request, StatusCode};
use koenig_damico_planner_server::app;
use serde_json::Value;
use tower::ServiceExt; // for `oneshot`

const WORKED_EXAMPLE_JSON: &str = r#"{
    "chief": {"a": 25000000.0, "e": 0.7, "i": 40.0, "raan": 358.0, "argp": 0.0, "mean_anom": 180.0},
    "t_i": 0.0,
    "t_f": 117990.0,
    "dt": 30.0,
    "w_metres": [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
    "cost": {"type": "piecewise"}
}"#;

/// Same valid window as the worked example but a tiny `dt` → ~1.18e8 grid points,
/// far past `MAX_GRID_POINTS`. The grid cap (audit B1) must reject this as a 400
/// at the HTTP boundary *before* any solve allocation.
const OVERSIZED_GRID_JSON: &str = r#"{
    "chief": {"a": 25000000.0, "e": 0.7, "i": 40.0, "raan": 358.0, "argp": 0.0, "mean_anom": 180.0},
    "t_i": 0.0,
    "t_f": 117990.0,
    "dt": 0.001,
    "w_metres": [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
    "cost": {"type": "norm2"}
}"#;

#[tokio::test]
async fn health_ok() {
    let req = Request::builder()
        .method("GET")
        .uri("/health")
        .body(Body::empty())
        .unwrap();

    let resp = app().oneshot(req).await.unwrap();
    assert_eq!(resp.status(), StatusCode::OK);

    let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
        .await
        .unwrap();
    let v: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(v["status"], "ok");
}

#[tokio::test]
async fn post_solve_golden_returns_200_within_bands() {
    let req = Request::builder()
        .method("POST")
        .uri("/solve")
        .header(header::CONTENT_TYPE, "application/json")
        .body(Body::from(WORKED_EXAMPLE_JSON))
        .unwrap();

    let resp = app().oneshot(req).await.unwrap();
    assert_eq!(resp.status(), StatusCode::OK);

    let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
        .await
        .unwrap();
    let v: Value = serde_json::from_slice(&bytes).unwrap();

    let maneuvers = v["maneuvers"].as_array().expect("maneuvers array");
    assert!(
        (1..=6).contains(&maneuvers.len()),
        "expected 1..=6 maneuvers, got {}",
        maneuvers.len()
    );

    let total_dv = v["total_dv"].as_f64().expect("total_dv number");
    assert!(
        (0.078..0.083).contains(&total_dv),
        "total_dv {total_dv} out of band (0.078, 0.083)"
    );

    let residual = v["residual"].as_f64().expect("residual number");
    assert!(residual < 1e-3, "residual {residual} must be < 1e-3");

    let lambda = v["lambda"].as_array().expect("lambda array");
    assert_eq!(lambda.len(), 6, "lambda must have 6 elements");
}

#[tokio::test]
async fn post_solve_malformed_returns_400_bad_request() {
    let req = Request::builder()
        .method("POST")
        .uri("/solve")
        .header(header::CONTENT_TYPE, "application/json")
        .body(Body::from("{ not json"))
        .unwrap();

    let resp = app().oneshot(req).await.unwrap();
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);

    let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
        .await
        .unwrap();
    let v: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(v["kind"], "bad_request");
}

#[tokio::test]
async fn post_solve_oversized_grid_returns_400_bad_request() {
    let req = Request::builder()
        .method("POST")
        .uri("/solve")
        .header(header::CONTENT_TYPE, "application/json")
        .body(Body::from(OVERSIZED_GRID_JSON))
        .unwrap();

    let resp = app().oneshot(req).await.unwrap();
    assert_eq!(resp.status(), StatusCode::BAD_REQUEST);

    let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
        .await
        .unwrap();
    let v: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(v["kind"], "bad_request");
    assert!(
        v["message"].as_str().unwrap().contains("grid"),
        "message should name the grid-size cap, got {:?}",
        v["message"]
    );
}
