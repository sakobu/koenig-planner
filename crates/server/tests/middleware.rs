//! Middleware behavior: permissive CORS header + request-body size cap.

use axum::body::Body;
use axum::http::{header, Request, StatusCode};
use koenig_damico_planner_server::app;
use tower::ServiceExt;

#[tokio::test]
async fn cors_allows_any_origin() {
    let req = Request::builder()
        .method("GET")
        .uri("/health")
        .header(header::ORIGIN, "https://example.com")
        .body(Body::empty())
        .unwrap();

    let resp = app().oneshot(req).await.unwrap();
    assert_eq!(resp.status(), StatusCode::OK);
    assert_eq!(
        resp.headers()
            .get(header::ACCESS_CONTROL_ALLOW_ORIGIN)
            .map(|v| v.to_str().unwrap().to_owned()),
        Some("*".to_owned())
    );
}

#[tokio::test]
async fn oversized_body_is_rejected_413() {
    // 70 KiB exceeds the 64 KiB cap; Content-Length short-circuits at the layer.
    let big = "x".repeat(70 * 1024);
    let req = Request::builder()
        .method("POST")
        .uri("/solve")
        .header(header::CONTENT_TYPE, "application/json")
        .body(Body::from(big))
        .unwrap();

    let resp = app().oneshot(req).await.unwrap();
    assert_eq!(resp.status(), StatusCode::PAYLOAD_TOO_LARGE);
}
