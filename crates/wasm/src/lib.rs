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
mod roe_track;

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

/// Plan a maneuver set from a typed request. The outcome (success or error) is
/// returned as a [`dto::SolveOutcome`] value so the error type is visible in
/// TypeScript.
///
/// # Never throws — for schema-conforming input
/// The body never panics or throws; every failure is an `Err` variant of the
/// returned value. The one throw is *outside* the body: called from untyped JS
/// with a value that is not a valid `SolveRequest` (missing / mistyped fields),
/// wasm-bindgen rejects it at the argument-deserialization boundary before this
/// function runs. A TypeScript caller cannot reach that path — the typed
/// `SolveRequest` parameter makes it a compile error.
///
/// # Unknown fields
/// The typed path (serde-wasm-bindgen) silently **ignores** unknown top-level or
/// `params` keys, whereas [`solve_json`] and the HTTP / Python frontends reject
/// them (`#[serde(deny_unknown_fields)]`). Adding `deny_unknown_fields` to the
/// wasm mirror DTOs is a no-op — serde-wasm-bindgen never surfaces unknown JS
/// properties — but a TypeScript `SolveRequest` forbids unknown keys at compile
/// time regardless.
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

/// Batch min-fuel dual over `base`'s window for many targets. Returns the gauge
/// `c*` (m/s; `None` if unreachable) and dual normal `λ` per target — the
/// reachable-set / Δv cost-map primitive. Never returns maneuvers.
#[wasm_bindgen]
pub fn sweep_dual(req: dto::SweepRequest) -> dto::SweepOutcome {
    let api_base: koenig_damico_planner_api::SolveRequest = (&req.base).into();
    match koenig_damico_planner_api::sweep(&api_base, &req.w_list) {
        Ok(points) => dto::SweepOutcome::Ok {
            value: points.into_iter().map(Into::into).collect(),
        },
        Err(e) => dto::SweepOutcome::Err { error: e.into() },
    }
}

/// Batch min-fuel **primal** solve over `base`'s window for many targets — the
/// reachable-set / Δv cost-map engine. Returns cost `c*`, dual `λ`, feasibility,
/// refine `iterations`, confidence `residual`, and burn count per target.
/// Reachability is the `feasible` field (not a `residual` threshold); on feasible
/// cells `residual` is ~machine-zero when the cell is clean and well-conditioned
/// but rises (to ~1e-7) as the extract/recovery step becomes numerically
/// ill-conditioned — a conditioning proxy, not a reachability metric. Never
/// returns maneuvers.
#[wasm_bindgen]
pub fn sweep_solve(req: dto::SweepRequest) -> dto::SweepSolveOutcome {
    let api_base: koenig_damico_planner_api::SolveRequest = (&req.base).into();
    match koenig_damico_planner_api::sweep_solve(&api_base, &req.w_list) {
        Ok(points) => dto::SweepSolveOutcome::Ok {
            value: points.into_iter().map(Into::into).collect(),
        },
        Err(e) => dto::SweepSolveOutcome::Err { error: e.into() },
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
