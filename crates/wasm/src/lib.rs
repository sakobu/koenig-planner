//! WebAssembly bindings for the Koenig-D'Amico maneuver planner.
//!
//! Exposes a typed [`solve`] (outcome modeled as a value so the error type
//! survives into the generated `.d.ts`) and a string [`solve_json`] escape
//! hatch. All wasm-bindgen / tsify concerns live in this leaf crate; the
//! shared `crates/api` facade stays platform-agnostic.

mod convert;
mod dto;
mod frames;
mod geometry;

pub use dto::*;

use wasm_bindgen::prelude::*;

/// Install a panic hook so Rust panics surface as readable console errors.
#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

/// Crate version (smoke test + demo footer).
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// Plan a maneuver set from a typed request. NEVER throws: the outcome
/// (success or error) is returned as a [`dto::SolveOutcome`] value so the
/// error type is visible in TypeScript.
#[wasm_bindgen]
pub fn solve(req: dto::SolveRequest) -> dto::SolveOutcome {
    let api_req: koenig_damico_planner_api::SolveRequest = (&req).into();
    match koenig_damico_planner_api::run(api_req) {
        Ok(resp) => match geometry::chief_geometry(&req, &resp) {
            Ok(geom) => dto::SolveOutcome::Ok {
                value: (resp, geom).into(),
            },
            // The solve succeeded; only presentation geometry failed. Deputy-
            // derived fields degrade in place (see geometry::chief_geometry), so
            // this arm is unreachable for a solved request — a residual error is
            // an internal invariant break, not a solver or input failure.
            Err(e) => dto::SolveOutcome::Err {
                error: dto::ApiError {
                    kind: dto::ApiErrorKind::Internal,
                    message: format!("presentation geometry: {e}"),
                },
            },
        },
        Err(e) => dto::SolveOutcome::Err { error: e.into() },
    }
}

/// String-in / string-out escape hatch (delegates to `api::run_json`). Returns
/// the response JSON on success; throws the serialized `ApiError` JSON on failure.
#[wasm_bindgen]
pub fn solve_json(json: &str) -> Result<String, JsValue> {
    koenig_damico_planner_api::run_json(json).map_err(|e| {
        JsValue::from_str(&serde_json::to_string(&e).unwrap_or_else(|_| e.to_string()))
    })
}
