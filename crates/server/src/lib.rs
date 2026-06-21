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
use koenig_damico_planner_api::{run, ApiError, ApiErrorKind, SolveRequest, SolveResponse};
use std::any::Any;
use std::time::Duration;
use tower::limit::ConcurrencyLimitLayer;
use tower::ServiceBuilder;
use tower_http::catch_panic::CatchPanicLayer;
use tower_http::cors::CorsLayer;
use tower_http::limit::RequestBodyLimitLayer;
use tower_http::timeout::TimeoutLayer;
use tower_http::trace::TraceLayer;

/// Default cap on simultaneous in-flight solves. Worst-case memory ≈
/// `DEFAULT_MAX_CONCURRENCY × (MAX_GRID_POINTS × 144 B)` ≈ 64 × 14 MB ≈ 900 MB.
const DEFAULT_MAX_CONCURRENCY: usize = 64;
/// Default per-request timeout (seconds). A solve is sub-second even at the grid
/// cap, so this only sheds genuinely-stuck requests.
const DEFAULT_TIMEOUT_SECS: u64 = 10;

const TIMEOUT_ENV: &str = "KOENIG_PLANNER_TIMEOUT_SECS";
const CONCURRENCY_ENV: &str = "KOENIG_PLANNER_MAX_CONCURRENCY";

/// Resolved transport-hardening limits applied to the router.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ServerConfig {
    pub(crate) max_concurrency: usize,
    pub(crate) timeout: Duration,
}

/// Pure config resolver: parse optional env-var strings, falling back to the
/// defaults on absent / unparseable / non-positive values. Never panics.
pub(crate) fn parse_config(
    timeout_secs: Option<String>,
    max_concurrency: Option<String>,
) -> ServerConfig {
    let timeout_secs = timeout_secs
        .and_then(|s| s.trim().parse::<u64>().ok())
        .filter(|&v| v > 0)
        .unwrap_or(DEFAULT_TIMEOUT_SECS);
    let max_concurrency = max_concurrency
        .and_then(|s| s.trim().parse::<usize>().ok())
        .filter(|&v| v > 0)
        .unwrap_or(DEFAULT_MAX_CONCURRENCY);
    ServerConfig {
        max_concurrency,
        timeout: Duration::from_secs(timeout_secs),
    }
}

/// Resolve [`ServerConfig`] from the process environment.
fn config_from_env() -> ServerConfig {
    parse_config(
        std::env::var(TIMEOUT_ENV).ok(),
        std::env::var(CONCURRENCY_ENV).ok(),
    )
}

/// Liveness handler: returns `200` with a small status body.
async fn health() -> Json<serde_json::Value> {
    Json(serde_json::json!({ "status": "ok" }))
}

/// Build the generic internal-fault [`ApiError`] for a caught panic, logging the
/// payload server-side. The client never receives the payload (no info leak).
fn internal_error(panic: &(dyn Any + Send)) -> ApiError {
    let detail = panic
        .downcast_ref::<&str>()
        .map(|s| (*s).to_owned())
        .or_else(|| panic.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "non-string panic payload".to_owned());
    tracing::error!(panic = %detail, "request handler panicked");
    ApiError {
        kind: ApiErrorKind::Internal,
        message: "internal server error".into(),
    }
}

/// `CatchPanicLayer` handler: convert a caught panic into the uniform 500.
fn handle_panic(err: Box<dyn Any + Send + 'static>) -> Response {
    AppError::from(internal_error(&*err)).into_response()
}

/// Map an [`ApiErrorKind`] to the response status code. Exhaustive: a new kind
/// is a compile error here, never a silent fall-through.
pub(crate) fn status_for_kind(kind: ApiErrorKind) -> StatusCode {
    match kind {
        ApiErrorKind::BadRequest => StatusCode::BAD_REQUEST,
        ApiErrorKind::Solver => StatusCode::UNPROCESSABLE_ENTITY,
        ApiErrorKind::Internal => StatusCode::INTERNAL_SERVER_ERROR,
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
                    kind: ApiErrorKind::BadRequest,
                    message: rejection.body_text(),
                },
            }),
        }
    }
}

/// Plan a maneuver set. Body is a `SolveRequest`; response is a `SolveResponse`.
///
/// The solve is synchronous CPU work, so it runs on the blocking pool: this keeps
/// the async reactor free for liveness (`/health`) and lets the `TimeoutLayer`
/// actually elapse (a non-yielding handler would never let the timer fire).
///
/// On timeout the layer drops this future — releasing the concurrency permit and
/// returning 408 — but the abandoned `spawn_blocking` task is uncancellable and
/// runs to completion on the blocking pool. That is safe only because the
/// `MAX_GRID_POINTS` cap (in `api::run`) bounds each solve to tens of ms / ~14 MB.
///
/// A solve-task panic surfaces as a `JoinError` → uniform 500 `{kind:"internal"}`.
async fn solve(ApiJson(req): ApiJson<SolveRequest>) -> Result<Json<SolveResponse>, AppError> {
    let resp = tokio::task::spawn_blocking(move || run(req))
        .await
        // spawn_blocking tasks cannot be cancelled, so a JoinError is always a
        // panic; route it through the same logged internal-error path. The outer
        // `?` maps ApiError→AppError (From); the inner `?` maps the api ApiError.
        .map_err(|join_err| internal_error(&*join_err.into_panic()))??;
    Ok(Json(resp))
}

/// The application routes (no middleware).
fn routes() -> Router {
    Router::new()
        .route("/health", get(health))
        .route("/solve", post(solve))
}

/// Apply transport-hardening middleware to a router.
///
/// Layer order (outermost → innermost): Trace, CORS, CatchPanic, ConcurrencyLimit,
/// Timeout, BodyLimit. CatchPanic sits *inside* CORS and Trace so a panic-derived
/// 500 still gets the CORS header and a trace log, and *outside* the handler and
/// the remaining middleware so it catches their panics. (Solve-closure panics are
/// handled separately, via the `spawn_blocking` JoinError arm.) The router's error
/// type stays `Infallible`, so no `HandleErrorLayer` is needed.
///
/// Note: `RequestBodyLimitLayer` is applied via `Router::layer` rather than inside
/// `ServiceBuilder` because `tower_http::limit::ResponseBody` does not implement
/// `Default`, which `ConcurrencyLimit`'s tower composition requires. The ordering
/// is preserved: BodyLimit is innermost (applied first by `Router::layer`), then
/// the `ServiceBuilder` wraps it with Timeout, ConcurrencyLimit, CatchPanic, CORS,
/// and Trace.
pub(crate) fn harden(router: Router, cfg: ServerConfig) -> Router {
    router.layer(RequestBodyLimitLayer::new(64 * 1024)).layer(
        ServiceBuilder::new()
            .layer(TraceLayer::new_for_http())
            .layer(CorsLayer::permissive())
            .layer(CatchPanicLayer::custom(handle_panic))
            .layer(ConcurrencyLimitLayer::new(cfg.max_concurrency))
            .layer(TimeoutLayer::with_status_code(
                StatusCode::REQUEST_TIMEOUT, // (status, duration)
                cfg.timeout,
            )),
    )
}

/// Build the application router with transport-hardening middleware.
pub fn app() -> Router {
    harden(routes(), config_from_env())
}

#[cfg(test)]
mod tests {
    use super::*;
    use axum::body::Body;
    use axum::extract::FromRequest;
    use axum::http::{header, Request, StatusCode};
    use axum::response::IntoResponse;
    use axum::routing::get;
    use koenig_damico_planner_api::ApiError;
    use serde_json::Value;
    use tower::ServiceExt;

    #[test]
    fn app_error_maps_bad_request_to_400() {
        let e: AppError = ApiError {
            kind: ApiErrorKind::BadRequest,
            message: "x".into(),
        }
        .into();
        assert_eq!(e.into_response().status(), StatusCode::BAD_REQUEST);
    }

    #[test]
    fn app_error_maps_solver_to_422() {
        let e: AppError = ApiError {
            kind: ApiErrorKind::Solver,
            message: "x".into(),
        }
        .into();
        assert_eq!(e.into_response().status(), StatusCode::UNPROCESSABLE_ENTITY);
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

    #[test]
    fn app_error_maps_internal_to_500() {
        let e: AppError = ApiError {
            kind: ApiErrorKind::Internal,
            message: "x".into(),
        }
        .into();
        assert_eq!(
            e.into_response().status(),
            StatusCode::INTERNAL_SERVER_ERROR
        );
    }

    #[test]
    fn parse_config_uses_defaults_when_absent_or_invalid() {
        // Absent → defaults.
        let c = parse_config(None, None);
        assert_eq!(c.max_concurrency, DEFAULT_MAX_CONCURRENCY);
        assert_eq!(
            c.timeout,
            std::time::Duration::from_secs(DEFAULT_TIMEOUT_SECS)
        );

        // Unparseable → defaults (must not panic).
        let c = parse_config(Some("abc".into()), Some("-5".into()));
        assert_eq!(c.max_concurrency, DEFAULT_MAX_CONCURRENCY);
        assert_eq!(
            c.timeout,
            std::time::Duration::from_secs(DEFAULT_TIMEOUT_SECS)
        );

        // Zero → defaults (0 concurrency would deadlock; 0 s would reject all).
        let c = parse_config(Some("0".into()), Some("0".into()));
        assert_eq!(c.max_concurrency, DEFAULT_MAX_CONCURRENCY);
        assert_eq!(
            c.timeout,
            std::time::Duration::from_secs(DEFAULT_TIMEOUT_SECS)
        );
    }

    #[test]
    fn parse_config_reads_valid_values() {
        let c = parse_config(Some("30".into()), Some("128".into()));
        assert_eq!(c.timeout, std::time::Duration::from_secs(30));
        assert_eq!(c.max_concurrency, 128);
    }

    #[tokio::test]
    async fn handle_panic_returns_uniform_internal_500() {
        let resp = handle_panic(Box::new("boom"));
        assert_eq!(resp.status(), StatusCode::INTERNAL_SERVER_ERROR);
        let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
            .await
            .unwrap();
        let v: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(v["kind"], "internal");
        assert_eq!(v["message"], "internal server error");
    }

    #[tokio::test]
    async fn caught_panic_maps_to_uniform_500_with_cors_header() {
        let cfg = parse_config(None, None);
        // An explicit `-> Response` return type keeps the diverging handler from
        // resolving to `!` (which fails `!: IntoResponse` under Rust 1.92's
        // `rust_2024_compatibility` never-type-fallback lint); the panic itself is
        // what the test exercises.
        async fn boom() -> Response {
            panic!("kaboom")
        }
        let panicking = Router::new().route("/boom", get(boom));
        let app = harden(panicking, cfg);

        let req = Request::builder()
            .method("GET")
            .uri("/boom")
            .header(header::ORIGIN, "https://example.com")
            .body(Body::empty())
            .unwrap();

        let resp = app.oneshot(req).await.unwrap();
        assert_eq!(resp.status(), StatusCode::INTERNAL_SERVER_ERROR);
        assert_eq!(
            resp.headers()
                .get(header::ACCESS_CONTROL_ALLOW_ORIGIN)
                .map(|v| v.to_str().unwrap().to_owned()),
            Some("*".to_owned()),
            "the synthesized 500 must still carry the permissive CORS header"
        );

        let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
            .await
            .unwrap();
        let v: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(v["kind"], "internal");
        assert_eq!(v["message"], "internal server error");
    }
}
