//! tsify-derived mirror DTOs for the wasm boundary.
//!
//! Field-for-field mirrors of the `crates/api` DTOs (request types add
//! `from_wasm_abi`, response types add `into_wasm_abi`), plus the presentation
//! `geometry` block and the `SolveOutcome` tagged-union result. The `From` impls
//! in `convert.rs` keep these in lock-step with `api` at compile time, on both
//! the request and response paths.

use serde::{Deserialize, Serialize};
use tsify::Tsify;

// ── Request side (deserialized from JS) ────────────────────────────────────

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(from_wasm_abi)]
pub struct OrbitDto {
    pub a: f64,
    pub e: f64,
    pub i: f64,
    pub raan: f64,
    pub argp: f64,
    pub mean_anom: f64,
}

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(from_wasm_abi)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum CostSpec {
    Norm2,
    FaceMax,
    Piecewise {
        #[tsify(optional)]
        period: Option<f64>,
        #[tsify(optional)]
        t_perigee0: Option<f64>,
    },
}

#[derive(Tsify, Serialize, Deserialize, Clone, Default)]
#[tsify(from_wasm_abi)]
pub struct SolveParamsDto {
    #[tsify(optional)]
    pub n_coarse: Option<usize>,
    #[tsify(optional)]
    pub n_init: Option<usize>,
    #[tsify(optional)]
    pub eps_cost: Option<f64>,
    #[tsify(optional)]
    pub eps_remove: Option<f64>,
}

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(from_wasm_abi)]
pub struct SolveRequest {
    pub chief: OrbitDto,
    pub t_i: f64,
    pub t_f: f64,
    pub dt: f64,
    pub w_metres: [f64; 6],
    pub cost: CostSpec,
    #[tsify(optional)]
    pub params: Option<SolveParamsDto>,
    #[tsify(optional)]
    pub initial_times: Option<Vec<f64>>,
}

// ── Response side (serialized to JS) ───────────────────────────────────────

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ManeuverDto {
    pub t: f64,
    pub dv: [f64; 3],
}

/// Presentation-only geometry derived from the chief via the core's verified
/// Kepler solver (see `geometry.rs`). `maneuver_nu[j]` is the true anomaly at
/// maneuver `j`; `perigee_window` (piecewise only) is the FaceMax band `[lo, hi]`
/// in true anomaly.
#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ChiefGeometry {
    pub a: f64,
    pub e: f64,
    pub maneuver_nu: Vec<f64>,
    #[tsify(optional)]
    pub perigee_window: Option<[f64; 2]>,
}

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct SolveResponse {
    pub maneuvers: Vec<ManeuverDto>,
    /// Total fuel cost [m/s]: the minimized objective (the paper's "delta-v
    /// cost" c*) — Σ‖Δv‖₂ under the L2 cost, the polytope gauge Σθ under FaceMax.
    pub total_dv: f64,
    pub iterations: usize,
    pub residual: f64,
    pub lambda: [f64; 6],
    pub geometry: ChiefGeometry,
}

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ApiError {
    pub kind: String,
    pub message: String,
}

/// Outcome modeled as a value so the error type survives into the `.d.ts`
/// (a wasm `Result` would erase `Err` into an untyped JS throw).
#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
#[serde(tag = "status", rename_all = "lowercase")]
pub enum SolveOutcome {
    Ok { value: SolveResponse },
    Err { error: ApiError },
}
