//! HTTP service for the Koenig-D'Amico maneuver planner.
//!
//! Thin wrapper over `koenig_damico_planner_api`: this crate owns transport and
//! runtime concerns only and never touches the generic solver.

use axum::extract::rejection::JsonRejection;
use axum::extract::{FromRequest, Request};
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
use koenig_damico_planner_api::{run, ApiError, SolveRequest, SolveResponse};
use tower_http::cors::CorsLayer;
use tower_http::limit::RequestBodyLimitLayer;
use tower_http::trace::TraceLayer;

/// Liveness handler: returns `200` with a small status body.
async fn health() -> Json<serde_json::Value> {
    Json(serde_json::json!({ "status": "ok" }))
}

/// Map an `ApiError` `kind` to the response status code.
pub(crate) fn status_for_kind(kind: &str) -> StatusCode {
    match kind {
        "bad_request" => StatusCode::BAD_REQUEST,
        "solver" => StatusCode::UNPROCESSABLE_ENTITY,
        _ => StatusCode::INTERNAL_SERVER_ERROR,
    }
}

/// Response wrapper that carries an explicit status plus a serializable
/// `ApiError` body, so every error path returns the same `{kind, message}` shape.
pub(crate) struct AppError {
    pub(crate) status: StatusCode,
    pub(crate) body: ApiError,
}

impl From<ApiError> for AppError {
    fn from(body: ApiError) -> Self {
        Self {
            status: status_for_kind(body.kind),
            body,
        }
    }
}

impl IntoResponse for AppError {
    fn into_response(self) -> Response {
        (self.status, Json(self.body)).into_response()
    }
}

/// Like `axum::Json`, but rejections are reshaped to the `{kind, message}`
/// contract (with axum's status preserved) instead of axum's default text body.
pub(crate) struct ApiJson<T>(pub T);

impl<S, T> FromRequest<S> for ApiJson<T>
where
    Json<T>: FromRequest<S, Rejection = JsonRejection>,
    S: Send + Sync,
{
    type Rejection = AppError;

    async fn from_request(req: Request, state: &S) -> Result<Self, Self::Rejection> {
        match Json::<T>::from_request(req, state).await {
            Ok(Json(value)) => Ok(ApiJson(value)),
            Err(rejection) => Err(AppError {
                status: rejection.status(),
                body: ApiError {
                    kind: "bad_request",
                    message: rejection.body_text(),
                },
            }),
        }
    }
}

/// Plan a maneuver set. Body is a `SolveRequest`; response is a `SolveResponse`.
async fn solve(ApiJson(req): ApiJson<SolveRequest>) -> Result<Json<SolveResponse>, AppError> {
    Ok(Json(run(req)?))
}

/// Build the application router.
pub fn app() -> Router {
    Router::new()
        .route("/health", get(health))
        .route("/solve", post(solve))
        .layer(TraceLayer::new_for_http())
        .layer(CorsLayer::permissive())
        .layer(RequestBodyLimitLayer::new(64 * 1024))
}

#[cfg(test)]
mod tests {
    use super::*;
    use axum::body::Body;
    use axum::extract::FromRequest;
    use axum::http::{header, Request, StatusCode};
    use axum::response::IntoResponse;
    use koenig_damico_planner_api::ApiError;
    use serde_json::Value;

    #[test]
    fn app_error_maps_bad_request_to_400() {
        let e: AppError = ApiError {
            kind: "bad_request",
            message: "x".into(),
        }
        .into();
        assert_eq!(e.into_response().status(), StatusCode::BAD_REQUEST);
    }

    #[test]
    fn app_error_maps_solver_to_422() {
        let e: AppError = ApiError {
            kind: "solver",
            message: "x".into(),
        }
        .into();
        assert_eq!(e.into_response().status(), StatusCode::UNPROCESSABLE_ENTITY);
    }

    #[test]
    fn app_error_maps_unknown_to_500() {
        let e: AppError = ApiError {
            kind: "mystery",
            message: "x".into(),
        }
        .into();
        assert_eq!(
            e.into_response().status(),
            StatusCode::INTERNAL_SERVER_ERROR
        );
    }

    #[tokio::test]
    async fn api_json_rejection_is_bad_request_kind_with_preserved_status() {
        // Malformed JSON syntax -> axum's JsonSyntaxError -> 400, body reshaped.
        let req = Request::builder()
            .method("POST")
            .header(header::CONTENT_TYPE, "application/json")
            .body(Body::from("{ not json"))
            .unwrap();

        let rejected = <ApiJson<Value> as FromRequest<()>>::from_request(req, &())
            .await
            .err()
            .expect("malformed json must be rejected");

        let resp = rejected.into_response();
        assert_eq!(resp.status(), StatusCode::BAD_REQUEST);

        let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
            .await
            .unwrap();
        let v: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(v["kind"], "bad_request");
        assert!(v["message"].is_string());
    }
}
