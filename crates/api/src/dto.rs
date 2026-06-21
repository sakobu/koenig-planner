//! Serde request/response DTOs for the planner wire contract.
//!
//! These mirror the core types but decouple the wire/JSON contract from the
//! nalgebra-based domain types: requests carry degrees/metres, responses carry
//! plain arrays. `convert.rs` holds the (field-exhaustive) conversions.

use serde::{Deserialize, Serialize};

// ── Request DTOs ────────────────────────────────────────────────────────────

/// Chief orbit definition.  Angles are in **degrees** (converted to radians
/// server-side); `a` is in **metres**.
///
/// These are the six mean Keplerian elements `[a, e, i, Ω, ω, M]` as used
/// throughout Koenig & D'Amico (2020).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitDto {
    /// Semimajor axis `[m]`.
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination `[deg]`.
    pub i: f64,
    /// Right ascension of the ascending node, Ω `[deg]`.
    pub raan: f64,
    /// Argument of perigee, ω `[deg]`.
    pub argp: f64,
    /// Mean anomaly, M `[deg]`.
    pub mean_anom: f64,
}

/// Which cost model to apply at each maneuver time.
///
/// For `Piecewise`, `period` defaults to the chief's Keplerian orbit period
/// (`2π / n`) when omitted — supplying a period unrelated to the chief
/// silently misaligns the perigee windows, so prefer the default.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum CostSpec {
    /// L2 norm (isotropic cost).
    Norm2,
    /// FaceMax gauge (fuel-optimal for an impulsive thruster set).
    FaceMax,
    /// Piecewise eq.-49 selector: FaceMax near perigee, Norm2 elsewhere.
    Piecewise {
        /// Orbit period `[s]`.  When `None`, derived as `2π / n` from the
        /// chief — strongly preferred so the perigee windows align correctly.
        #[serde(default)]
        period: Option<f64>,
        /// First perigee-passage epoch `[s]`.  When `None`, defaults to
        /// `period / 2` (apogee-at-`t = 0` assumption, matching the worked
        /// example where `M₀ = 180°`).
        #[serde(default)]
        t_perigee0: Option<f64>,
    },
}

/// Solver tuning knobs.  Every field is optional; missing fields fall back to
/// [`SolveParams::default`](crate::core::SolveParams::default) (Table III of Koenig & D'Amico 2020).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolveParamsDto {
    /// Coarse-sample count `|T^d|` for Algorithm 1 (default 20).
    #[serde(default)]
    pub n_coarse: Option<usize>,
    /// Initial candidate-time count `n_init` (default 6).
    #[serde(default)]
    pub n_init: Option<usize>,
    /// Convergence tolerance `eps_cost` (default 0.01).
    #[serde(default)]
    pub eps_cost: Option<f64>,
    /// Slack-removal tolerance `eps_remove` (default 0.01).
    #[serde(default)]
    pub eps_remove: Option<f64>,
}

/// A full planning request.
///
/// Angles in [`OrbitDto`] are **degrees**; `w_metres` is in **metres**;
/// times are in **seconds**.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolveRequest {
    /// Chief mean absolute orbit (angles in degrees).
    pub chief: OrbitDto,
    /// Initial time `t_i` `[s]`.
    pub t_i: f64,
    /// Final time `t_f` `[s]`.
    pub t_f: f64,
    /// Grid step `dt` `[s]`.
    pub dt: f64,
    /// Target pseudostate in **metres**.  The server divides each component
    /// by `chief.a` to produce the dimensionless `w` passed to the planner,
    /// matching the nondimensionalisation in the worked example.
    pub w_metres: [f64; 6],
    /// Cost model selection.
    pub cost: CostSpec,
    /// Optional solver tuning (default = Table III).
    #[serde(default)]
    pub params: Option<SolveParamsDto>,
    /// Optional explicit initial candidate times for Algorithm 2 (bypasses
    /// Algorithm 1 when provided, enabling the paper's initialization study).
    #[serde(default)]
    pub initial_times: Option<Vec<f64>>,
}

// ── Response DTOs ───────────────────────────────────────────────────────────

/// A single impulsive maneuver in the RTN frame.
#[derive(Debug, Clone, Serialize)]
pub struct ManeuverDto {
    /// Application time `[s]`, measured from `t_i`.
    pub t: f64,
    /// Delta-v `[m/s]`, RTN components `[R, T, N]`.
    pub dv: [f64; 3],
}

/// Successful planner output.
#[derive(Debug, Clone, Serialize)]
pub struct SolveResponse {
    /// Ordered list of maneuvers.
    pub maneuvers: Vec<ManeuverDto>,
    /// Total fuel cost `[m/s]`: the minimized objective `Σⱼ f_{tⱼ}(Δvⱼ)` (the
    /// paper's "delta-v cost" `c*`). `Σ‖Δvⱼ‖₂` under the L2 cost; the polytope
    /// gauge `Σθ` (`≥` the net-Δv L2 norm) under FaceMax.
    pub total_dv: f64,
    /// Algorithm 2 iteration count.
    pub iterations: usize,
    /// Relative residual `‖w_err‖ / ‖w‖` of the pre-prune min-fuel solution.
    pub residual: f64,
    /// Optimal dual variable `λ_opt ∈ ℝ⁶`.
    pub lambda: [f64; 6],
}

// ── Error ───────────────────────────────────────────────────────────────────

/// Owned error that decouples the wire contract from [`PlannerError`](crate::core::PlannerError).
///
/// `kind` is an [`ApiErrorKind`]; serializes to `{"kind": …, "message": …}` so the
/// WASM/HTTP frontends can return it directly as a JSON error body.
#[derive(Debug, thiserror::Error, Serialize)]
#[error("{kind}: {message}")]
pub struct ApiError {
    /// Status class for HTTP/Python/WASM frontends.
    pub kind: ApiErrorKind,
    /// Human-readable description of what went wrong.
    pub message: String,
}

/// Status class for an [`ApiError`].
///
/// Serializes to a stable wire string via an **explicit per-variant**
/// `#[serde(rename)]` (deliberately *not* `rename_all`: the sibling cost/outcome
/// enums need single-word lowercase tags like `facemax`, whereas these need the
/// snake-cased `bad_request`; an explicit rename keeps each wire string literal
/// and local, immune to a `rename_all` change). `as_str` is the single source of
/// truth shared by [`Display`] and verified against serde by a wire-stability test.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ApiErrorKind {
    /// Invalid input / malformed request — the caller should fix the request.
    #[serde(rename = "bad_request")]
    BadRequest,
    /// Well-formed input but numerically unsolvable / solver failure.
    #[serde(rename = "solver")]
    Solver,
    /// Unexpected internal fault (e.g. a panic caught by the HTTP layer).
    #[serde(rename = "internal")]
    Internal,
}

impl ApiErrorKind {
    /// The stable wire string (`"bad_request"` / `"solver"` / `"internal"`).
    pub fn as_str(self) -> &'static str {
        match self {
            ApiErrorKind::BadRequest => "bad_request",
            ApiErrorKind::Solver => "solver",
            ApiErrorKind::Internal => "internal",
        }
    }
}

impl std::fmt::Display for ApiErrorKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}
