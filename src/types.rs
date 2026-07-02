//! Core value types: dimensions, pseudostate/maneuver, the uniform time grid,
//! solver parameters, the solver result, the conic-row placeholder, and errors.

use nalgebra::{SMatrix, SVector};
use thiserror::Error;

/// State dimension: 6 quasi-nonsingular relative orbital elements (ROEs).
pub const N: usize = 6;

/// Control dimension: 3 RTN Delta-v components (R, T, N).
pub const M: usize = 3;

/// A pseudostate / ROE vector in R^6 (dimensionless unless scaled by `a_c`).
pub type Pseudostate = SVector<f64, N>;

/// The dual variable lambda in R^6 (outward reachable-set normal).
pub type Dual = SVector<f64, N>;

/// An impulsive maneuver: a Delta-v `[m/s]` in the RTN frame applied at time `t` `[s]`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct Maneuver {
    /// Application time `[s]`: an absolute grid time `grid.time(k) = t_i + k·dt`
    /// (the same axis as the grid, equal to `t_i` at the first sample).
    pub t: f64,
    /// Delta-v [m/s], RTN components (R, T, N).
    pub dv: SVector<f64, M>,
}

/// A uniform time grid of samples `t_i + k·dt` lying within `[t_i, t_f]`, `dt > 0`.
///
/// The first sample is `t_i`; the last is the largest `t_i + k·dt` that does not
/// exceed `t_f`. `t_f` is itself a grid point exactly when the window length is a
/// whole multiple of `dt` (as in the paper's grids below); otherwise the final
/// sample falls short of `t_f` by less than `dt`. No candidate ever lands past
/// `t_f`, so the grid always respects the admissible domain (\[KD20\] eq. 5:
/// `T ⊆ [t_i, t_f]`).
///
/// The worked example (Table III) uses a 30 s grid over `[0, 117990]` -> 3934
/// candidate times; the Hunter cross-check uses a 10 s grid over `[0, 39000]`
/// -> 3901 candidate times.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TimeGrid {
    /// Initial time `t_i` `[s]`.
    pub t_i: f64,
    /// Final time `t_f` `[s]`.
    pub t_f: f64,
    /// Grid step `dt` `[s]`.
    pub dt: f64,
}

impl TimeGrid {
    /// Build a uniform grid with step `dt` over `[t_i, t_f]`.
    ///
    /// Ref: \[KD20\] worked-example control domain (uniform discretization of
    /// `[t_i, t_f]`), p. 10.
    ///
    /// # Errors
    /// Returns [`PlannerError::InvalidInput`] unless `t_i`, `t_f`, `dt` are all
    /// finite, `dt > 0`, and `t_f > t_i`. This is the only validating entry
    /// point: constructing `TimeGrid { .. }` via the public fields bypasses it,
    /// in which case `len`/`time`/`times` assume that same invariant.
    pub fn uniform(t_i: f64, t_f: f64, dt: f64) -> Result<Self, PlannerError> {
        if !t_i.is_finite() || !t_f.is_finite() || !dt.is_finite() || dt <= 0.0 || t_f <= t_i {
            return Err(PlannerError::InvalidInput(InvalidInputKind::Grid {
                t_i,
                t_f,
                dt,
            }));
        }
        Ok(Self { t_i, t_f, dt })
    }

    /// Number of grid points: the count of samples `t_i + k·dt` within
    /// `[t_i, t_f]` (always `>= 1`; includes `t_i`, and includes `t_f` only when
    /// the window length is a whole multiple of `dt`).
    ///
    /// Assumes the [`uniform`](Self::uniform) invariant (`dt > 0`, `t_f > t_i`,
    /// finite); on a hand-built `TimeGrid` violating it the `f64 -> usize` cast
    /// saturates.
    pub fn len(&self) -> usize {
        // floor (not round) so the last sample `t_i + (len-1)·dt` never lands past
        // `t_f` (\[KD20\] eq. 5: `T ⊆ [t_i, t_f]`). The relative tolerance restores
        // the exact endpoint on a commensurate window whose f64 division lands a
        // few ULP short of the whole ratio (e.g. `dt = (t_f - t_i)/n`), while
        // staying far below 1 so a genuinely short ratio is never pulled up.
        let ratio = (self.t_f - self.t_i) / self.dt;
        let tol = 1e-9 * ratio.abs().max(1.0);
        (ratio + tol).floor() as usize + 1
    }

    /// A grid always has at least one point; provided for lint-completeness.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The time `[s]` at grid index `idx`.
    pub fn time(&self, idx: usize) -> f64 {
        self.t_i + (idx as f64) * self.dt
    }

    /// Iterator over all grid times `[s]`.
    pub fn times(&self) -> impl Iterator<Item = f64> + '_ {
        (0..self.len()).map(move |i| self.time(i))
    }
}

/// Tunable parameters for the three-step algorithm (\[KD20\] p. 10 prose defaults).
///
/// The paper's Table III is the chief-orbit / pseudostate table; these solver
/// parameters come from the p. 10 prose ("20 times evenly distributed", "six
/// times", tolerance "selected as 0.01").
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SolveParams {
    /// Coarse-sample count `|T^d|` for Algorithm 1 initialization (p. 10 prose: 20).
    pub n_coarse: usize,
    /// Initial candidate-time count `n_init` (p. 10 prose: 6).
    pub n_init: usize,
    /// Convergence tolerance `eps_cost` (p. 10 prose: 0.01).
    pub eps_cost: f64,
    /// Slack-removal tolerance `eps_remove` (p. 10 prose: 0.01).
    pub eps_remove: f64,
}

impl Default for SolveParams {
    // Ref: [KD20] default params (T^d=20, T^est=6, eps=0.01), p. 10 prose.
    fn default() -> Self {
        Self {
            n_coarse: 20,
            n_init: 6,
            eps_cost: 0.01,
            eps_remove: 0.01,
        }
    }
}

/// The planner output: the maneuver set plus convergence diagnostics.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct Solution {
    /// Maneuvers `{t, Delta-v}`, one per optimal time in `T^opt`.
    pub maneuvers: Vec<Maneuver>,
    /// Total fuel cost `[m/s]`: the minimized objective `Σ_j f_{t_j}(Delta-v_j)`
    /// (the paper's "delta-v cost" `c*`; eq. 4). For the L2 cost this is
    /// `Σ ||Delta-v_j||`; for the FaceMax polytopic cost it is the gauge `Σ theta`
    /// (the sum of the tetrahedral-thruster firings), which is `>=` the L2 norm of
    /// the net Delta-v whenever a burn combines two or more vertices. This is the
    /// cost that was actually minimized — not `Σ ||Delta-v_j||` under FaceMax.
    ///
    /// Measured on the **full, pre-prune** solution (consistent with `residual`).
    pub total_dv: f64,
    /// Algorithm 2 iteration count.
    pub iterations: usize,
    /// Relative residual `||w_err|| / ||w||` of the **full, pre-prune** min-fuel
    /// solution over `T^opt` — the true reachability metric (approximately 0
    /// when `w` is reachable).
    ///
    /// Measured before interior-point dust is pruned from `maneuvers` (maneuvers
    /// below `1e-3` of the largest are dropped). Recomputing the residual from
    /// the returned, pruned `maneuvers` can therefore give a slightly larger
    /// value; it is bounded by the pruned mass and stays small. Use this field
    /// for the reachability check.
    pub residual: f64,
    /// Optimal dual lambda_opt.
    pub lambda: Dual,
}

/// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1` for one candidate time.
///
/// Linear rows encode `a^T lambda <= b`; SOC rows encode `||G lambda||_2 <= h`.
/// [`crate::solver::refine_socp()`] assembles these into clarabel cones.
#[derive(Debug, Clone, Default)]
pub struct ConicRows {
    /// Linear rows `(a, b)` with `a^T lambda <= b`.
    pub linear: Vec<(SVector<f64, N>, f64)>,
    /// Second-order-cone rows `(G, h)` with `||G lambda||_2 <= h`.
    pub soc: Vec<(SMatrix<f64, M, N>, f64)>,
}

/// Primal fuel generator for one maneuver in the direct min-fuel SOCP
/// (Algorithm 3). Describes how a Δv at one candidate time is built
/// from solver variables and how it is charged, mirroring the cost's unit
/// sublevel set:
///
/// * `Norm` — a free vector `v ∈ ℝᴹ` charged its L2 norm `‖v‖₂` (an `‖v‖ ≤ c`
///   second-order cone). This is the L2 cost.
/// * `Polytope(dirs)` — a nonnegative combination `Δv = Σₖ θₖ·dirs[k]`,
///   `θₖ ≥ 0`, charged `Σₖ θₖ` (a nonnegative-cone LP). The unit ball is
///   `conv{0, dirs…}`, so this is the gauge of that polytope; `dirs` are the
///   FaceMax `V_vertex` columns.
#[derive(Debug, Clone, PartialEq)]
pub enum FuelGenerator {
    /// L2: free `v ∈ ℝᴹ` charged `‖v‖₂`.
    Norm,
    /// Polytopic: `Δv = Σₖ θₖ·dirs[k]`, `θ ≥ 0`, charged `Σₖ θₖ`.
    Polytope(Vec<SVector<f64, M>>),
}

/// Classifies a [`PlannerError::InvalidInput`] and carries the offending value(s), so callers can
/// branch on the failure mode and read the diagnostic numbers without parsing a message string.
///
/// This enum is `#[non_exhaustive]`: match it with a trailing `_` arm so future kinds do not break
/// downstream code.
#[derive(Debug, Clone, PartialEq, Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[non_exhaustive]
pub enum InvalidInputKind {
    /// Time-grid parameters out of range.
    #[error("grid requires finite t_i, t_f, dt with dt > 0 and t_f > t_i (got t_i={t_i}, t_f={t_f}, dt={dt})")]
    Grid { t_i: f64, t_f: f64, dt: f64 },
    /// `n_init` or `n_coarse` was zero.
    #[error("n_init and n_coarse must both be >= 1 (got n_init={n_init}, n_coarse={n_coarse})")]
    SolverParams { n_init: usize, n_coarse: usize },
    /// Target pseudostate `w` was zero or non-finite.
    #[error("target pseudostate w must be nonzero and finite")]
    Target,
    /// Eccentricity outside the elliptic range.
    #[error("eccentricity must satisfy 0 <= e < 1 (elliptic), got e = {e}")]
    Eccentricity { e: f64 },
    /// Chief semimajor axis non-finite or non-positive.
    #[error("J2Roe: chief semimajor axis must be finite and positive (a > 0), got a = {a}")]
    ChiefSemimajorAxis { a: f64 },
    /// Chief inclination too close to 0 or pi.
    #[error("J2Roe: chief inclination must be bounded away from 0 and pi (B(t) has a 1/tan(i) singularity), got i = {i} rad")]
    ChiefInclination { i: f64 },
    /// Chief propagation window non-finite or `t_f <= t_i`.
    #[error("J2Roe: window must satisfy finite t_i, t_f and t_f > t_i (got t_i={t_i}, t_f={t_f})")]
    Window { t_i: f64, t_f: f64 },
    /// Orbital period non-finite/`<= 0`, or non-finite perigee epoch.
    #[error("Piecewise requires a finite period > 0 and a finite perigee epoch (got period={period}, t_perigee0={t_perigee0})")]
    Period { period: f64, t_perigee0: f64 },
    /// A maneuver budget was negative or NaN.
    #[error("extract_qp: budget must be non-negative, got {budget}")]
    Budget { budget: f64 },
    /// A candidate-time / generator set was empty or malformed.
    #[error("candidate-time set is empty or malformed")]
    EmptyCandidateSet,
    /// No maneuver directions supplied to the QP extraction.
    #[error("extract_qp: no maneuver directions")]
    NoDirections,
    /// No supplied initial time was finite and in range after snapping.
    #[error("solve_from_initial_times: no finite initial times in range")]
    NoInitialTimesInRange,
    /// An error fitting none of the planner's own preconditions — e.g. from an external
    /// `Dynamics`/`CostModel` implementation. Mirrors [`std::io::ErrorKind::Other`].
    #[error("{message}")]
    Other { message: String },
}

/// Errors surfaced by the planner.
///
/// This enum is `#[non_exhaustive]`: match it with a trailing `_` arm, or
/// classify it with [`PlannerError::class`], so future error categories do not
/// break downstream code.
#[derive(Debug, Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[non_exhaustive]
pub enum PlannerError {
    /// The convex (SOCP/QP) solver failed.
    #[error("convex solver failed: {0}")]
    SolverFailed(String),
    /// Algorithm 2 did not reach `max_t g <= 1 + eps_cost` within the iteration cap.
    #[error("refinement did not converge in {max_iters} iterations (max_t g = {achieved}, target <= {target})")]
    NotConverged {
        /// Iteration cap that was hit.
        max_iters: usize,
        /// Achieved `max_t g`.
        achieved: f64,
        /// Target threshold `1 + eps_cost`.
        target: f64,
    },
    /// The Kepler Newton iteration failed to converge.
    #[error("Kepler solve diverged for M = {m} rad, e = {e}")]
    KeplerDivergence {
        /// Mean anomaly `[rad]`.
        m: f64,
        /// Eccentricity.
        e: f64,
    },
    /// An input precondition was violated. The wrapped [`InvalidInputKind`] classifies the cause and
    /// carries the offending value(s); match on it to branch on the failure mode. Treat any
    /// `InvalidInput` uniformly as a "bad request — correct the inputs" signal.
    #[error("invalid input: {0}")]
    InvalidInput(InvalidInputKind),
}

/// Coarse, transport-agnostic category of a [`PlannerError`].
///
/// Lets a transport layer (HTTP status, process exit code, …) map an error
/// without matching every [`PlannerError`] variant. `#[non_exhaustive]`: match it
/// with a trailing `_` arm so a future category does not break downstream code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ErrorClass {
    /// Caller-fixable: the request was malformed; correct the inputs.
    /// (Wraps [`PlannerError::InvalidInput`].)
    InvalidInput,
    /// The request was well-formed but could not be solved — a numeric solver
    /// failure, non-convergence, or Kepler divergence.
    Unsolvable,
}

impl PlannerError {
    /// Classify this error into a coarse [`ErrorClass`] for transport mapping.
    ///
    /// The match is exhaustive (no wildcard), so a future `PlannerError` variant
    /// must be assigned a class here, inside this crate, at compile time. A
    /// downstream crate that maps `ErrorClass` instead keeps compiling and
    /// stays correct when a new `PlannerError` variant is added (it lands in an
    /// existing class) — which is what makes the `#[non_exhaustive]` on
    /// [`PlannerError`] non-breaking in practice.
    #[must_use]
    pub fn class(&self) -> ErrorClass {
        match self {
            PlannerError::InvalidInput(_) => ErrorClass::InvalidInput,
            PlannerError::SolverFailed(_)
            | PlannerError::NotConverged { .. }
            | PlannerError::KeplerDivergence { .. } => ErrorClass::Unsolvable,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn error_class_partitions_every_variant() {
        assert_eq!(
            PlannerError::InvalidInput(InvalidInputKind::Target).class(),
            ErrorClass::InvalidInput
        );
        assert_eq!(
            PlannerError::SolverFailed("x".into()).class(),
            ErrorClass::Unsolvable
        );
        assert_eq!(
            PlannerError::NotConverged {
                max_iters: 1,
                achieved: 2.0,
                target: 1.0,
            }
            .class(),
            ErrorClass::Unsolvable
        );
        assert_eq!(
            PlannerError::KeplerDivergence { m: 0.0, e: 0.7 }.class(),
            ErrorClass::Unsolvable
        );
    }

    // Ref: [KD20] eq. 51 (ROE state x = [da, dlambda, dex, dey, dix, diy]).
    #[test]
    fn dimensions_match_roe_and_rtn() {
        assert_eq!(N, 6);
        assert_eq!(M, 3);
    }

    // Ref: [KD20] worked-example control domain (30 s grid, 3934 times, t_f=117990 s),
    // p. 10.
    #[test]
    fn worked_example_grid_has_3934_times() {
        let g = TimeGrid::uniform(0.0, 117990.0, 30.0).unwrap();
        assert_eq!(g.len(), 3934);
        assert_abs_diff_eq!(g.time(0), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(g.time(g.len() - 1), 117990.0, epsilon = 1e-6);
        assert_eq!(g.times().count(), 3934);
    }

    // Ref: [H25] eq. 69 (10 s grid over [0, 39000]: 3901 times).
    #[test]
    fn hunter_grid_has_3901_times() {
        let g = TimeGrid::uniform(0.0, 39000.0, 10.0).unwrap();
        assert_eq!(g.len(), 3901);
    }

    // Ref: [KD20] eq. 5 (admissible domain T ⊆ [t_i, t_f]); eq. 11 (only impulses
    // with t_j ≤ t contribute to x(t)). On a window whose length is not a whole
    // multiple of `dt`, no grid time may land past `t_f` — such a candidate is
    // inadmissible and would be evaluated with a backward-extrapolated STM.
    #[test]
    fn grid_last_time_never_exceeds_t_f() {
        // Non-commensurate windows where round-based `len` used to overshoot.
        for (t_i, t_f, dt) in [
            (0.0, 117_990.0, 100.0), // ratio 1179.9 -> last was 118000 (+10 s)
            (0.0, 117_990.0, 29.0),  //             -> last was 118001 (+11 s)
            (0.0, 100.0, 40.0),      // ratio 2.5 (half away from zero) -> last 120
            (0.0, 5.0, 10.0),        // window shorter than dt -> only t_i fits
        ] {
            let g = TimeGrid::uniform(t_i, t_f, dt).unwrap();
            let last = g.time(g.len() - 1);
            assert!(
                last <= t_f,
                "uniform({t_i}, {t_f}, {dt}): last grid time {last} exceeds t_f {t_f}"
            );
            // Maximal count: the next sample must fall outside the window, so
            // flooring did not silently drop a point that still fits in [t_i, t_f].
            let beyond = g.time(g.len());
            assert!(
                beyond > t_f,
                "uniform({t_i}, {t_f}, {dt}): dropped an in-window sample at {beyond}"
            );
        }
    }

    // Ref: [KD20] default params (20, 6, 0.01, 0.01), p. 10 prose.
    #[test]
    fn default_params_match_paper() {
        let p = SolveParams::default();
        assert_eq!(p.n_coarse, 20);
        assert_eq!(p.n_init, 6);
        assert_abs_diff_eq!(p.eps_cost, 0.01, epsilon = 1e-12);
        assert_abs_diff_eq!(p.eps_remove, 0.01, epsilon = 1e-12);
    }

    #[test]
    fn conic_rows_default_is_empty() {
        let c = ConicRows::default();
        assert!(c.linear.is_empty());
        assert!(c.soc.is_empty());
    }

    #[test]
    fn uniform_rejects_degenerate_grids() {
        assert!(TimeGrid::uniform(0.0, 60.0, 0.0).is_err()); // dt = 0
        assert!(TimeGrid::uniform(0.0, 60.0, -1.0).is_err()); // dt < 0
        assert!(TimeGrid::uniform(0.0, 60.0, f64::NAN).is_err());
        assert!(TimeGrid::uniform(0.0, 60.0, f64::INFINITY).is_err());
        assert!(TimeGrid::uniform(60.0, 0.0, 1.0).is_err()); // t_f < t_i
        assert!(TimeGrid::uniform(0.0, 0.0, 1.0).is_err()); // zero-length window rejected (matches validate_inputs / J2Roe::new)
    }
}
