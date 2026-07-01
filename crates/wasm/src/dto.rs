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
#[serde(tag = "type")]
pub enum CostSpec {
    #[serde(rename = "norm2")]
    Norm2,
    #[serde(rename = "facemax")]
    FaceMax,
    #[serde(rename = "piecewise")]
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
    /// Coarse-sample count for Algorithm 1. Ignored when `initial_times` is
    /// supplied, since that path bypasses Algorithm 1.
    #[tsify(optional)]
    pub n_coarse: Option<usize>,
    /// Initial candidate-time count. Ignored when `initial_times` is supplied,
    /// since that path bypasses Algorithm 1.
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
    pub w_meters: [f64; 6],
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

/// A maneuver expressed in ECI for the 3D scene: the burn position on the chief
/// orbit and the executed Δv direction (RTN→ECI rotated). Presentation-only.
#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ManeuverEciDto {
    /// Burn position in ECI `[m]`.
    pub position_eci: [f64; 3],
    /// Executed Δv in ECI `[m/s]` (the thrust direction — distinct from the
    /// primer vector except under the `norm2` cost).
    pub dv_eci: [f64; 3],
}

/// A maneuver expressed in the chief RTN (relative-motion) frame for the 3D
/// scene: the deputy's relative burn position and the executed Δv. The RTN
/// analog of [`ManeuverEciDto`]. Presentation-only.
#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ManeuverRtnDto {
    /// Deputy relative burn position in the chief RTN frame `[m]`. A schematic
    /// anchor on the target-ROE deputy track (see `geometry.rs`); only the Δv
    /// direction is exact.
    pub position_rtn: [f64; 3],
    /// Executed Δv in the chief RTN frame `[m/s]` `[R, T, N]` — the native
    /// solver frame, echoed with no rotation.
    pub dv_rtn: [f64; 3],
}

/// Presentation-only geometry derived from the chief via the core's verified
/// Kepler solver (see `geometry.rs`). `maneuver_nu[j]` is the true anomaly at
/// maneuver `j`; `perigee_window` (piecewise only) is the FaceMax band `[lo, hi]`
/// in true anomaly. The `*_eci` fields are metric ECI samples for the 3D scene;
/// `relative_trajectory_rtn` is the deputy's relative orbit (meters) in the
/// chief RTN frame; `target_roe` echoes the request `w_meters`.
#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ChiefGeometry {
    pub a: f64,
    pub e: f64,
    pub maneuver_nu: Vec<f64>,
    #[tsify(optional)]
    pub perigee_window: Option<[f64; 2]>,
    /// Closed-loop chief-orbit samples in ECI `[m]` (orbit-shape curve).
    pub orbit_eci: Vec<[f64; 3]>,
    /// Chief position in ECI `[m]` at each `primer_times` sample (playback track).
    pub chief_track_eci: Vec<[f64; 3]>,
    /// Burn position + Δv direction in ECI, one per maneuver.
    pub maneuver_eci: Vec<ManeuverEciDto>,
    /// Burn position + native-RTN Δv per maneuver, in the chief RTN frame — the
    /// RTN analog of `maneuver_eci` for the relative-motion scene.
    pub maneuver_rtn: Vec<ManeuverRtnDto>,
    /// Primer vector in ECI at each `primer_times` sample (dimensionless dir).
    pub primer_eci: Vec<[f64; 3]>,
    /// Primer vector in the chief RTN frame at each `primer_times` sample — a
    /// presentation copy of the response `primer_rtn` (RTN analog of
    /// `primer_eci`), so the RTN scene consumes only this geometry block.
    pub primer_rtn: Vec<[f64; 3]>,
    /// ECI samples of the FaceMax perigee-window arc (piecewise cost only).
    #[tsify(optional)]
    pub perigee_arc_eci: Option<Vec<[f64; 3]>>,
    /// Deputy relative orbit in the chief RTN frame `[m]` (driven by `target_roe`).
    pub relative_trajectory_rtn: Vec<[f64; 3]>,
    /// Deputy position in the chief RTN frame `[m]` at each `primer_times` sample
    /// (the playback grid) — the deputy glyph that tracks the scrubber.
    pub deputy_track_rtn: Vec<[f64; 3]>,
    /// Echo of the request `w_meters` `[m]` = `[δa, δλ, δeₓ, δe_y, δiₓ, δi_y]·a`.
    pub target_roe: [f64; 6],
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
    /// Primer-vector history (the paper's Fig. 7), sampled at every grid time.
    /// The three arrays are parallel and equal length. `primer_times[k]` is the
    /// sample time `[s]`; `primer_magnitude[k]` is the dimensionless dual-gauge
    /// magnitude (`≤ 1 + eps_cost`, `≈ 1` at maneuver times); `primer_rtn[k]` is
    /// the primer vector `p(t) = Γᵀ(t)·λ` in RTN (not the executed thrust
    /// direction). Touch-1-away-from-a-burn ⇒ plan flexibility.
    pub primer_times: Vec<f64>,
    pub primer_magnitude: Vec<f64>,
    pub primer_rtn: Vec<[f64; 3]>,
    pub geometry: ChiefGeometry,
}

/// Status class for an [`ApiError`]. Mirrors `api::ApiErrorKind`; the explicit
/// per-variant `#[serde(rename)]` pins the wire tags (`bad_request`/`solver`/
/// `internal`) and types `ApiError.kind` in the generated `.d.ts` as the union
/// `"bad_request" | "solver" | "internal"` instead of a bare `string`.
#[derive(Tsify, Serialize, Deserialize, Clone, Copy, Debug, PartialEq, Eq)]
#[tsify(into_wasm_abi)]
pub enum ApiErrorKind {
    #[serde(rename = "bad_request")]
    BadRequest,
    #[serde(rename = "solver")]
    Solver,
    #[serde(rename = "internal")]
    Internal,
}

#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
pub struct ApiError {
    pub kind: ApiErrorKind,
    pub message: String,
}

/// Outcome modeled as a value so the error type survives into the `.d.ts`
/// (a wasm `Result` would erase `Err` into an untyped JS throw).
// SolveResponse is a wasm boundary type serialised immediately; stack size is
// not a concern here — suppress the large_enum_variant lint.
#[allow(clippy::large_enum_variant)]
#[derive(Tsify, Serialize, Deserialize, Clone)]
#[tsify(into_wasm_abi)]
#[serde(tag = "status")]
pub enum SolveOutcome {
    #[serde(rename = "ok")]
    Ok { value: SolveResponse },
    #[serde(rename = "err")]
    Err { error: ApiError },
}
