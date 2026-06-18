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

/// An impulsive maneuver: a Delta-v [m/s] in the RTN frame applied at time `t` [s].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Maneuver {
    /// Application time [s], measured from `t_i`.
    pub t: f64,
    /// Delta-v [m/s], RTN components (R, T, N).
    pub dv: SVector<f64, M>,
}

/// A uniform, endpoint-inclusive time grid over `[t_i, t_f]` with step `dt`.
///
/// The worked example (Table III) uses a 30 s grid over `[0, 117990]` -> 3934
/// candidate times; the Hunter cross-check uses a 10 s grid over `[0, 39000]`
/// -> 3901 candidate times.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeGrid {
    /// Initial time `t_i` [s].
    pub t_i: f64,
    /// Final time `t_f` [s].
    pub t_f: f64,
    /// Grid step `dt` [s].
    pub dt: f64,
}

impl TimeGrid {
    /// Build a uniform grid with step `dt` over `[t_i, t_f]`.
    pub fn uniform(t_i: f64, t_f: f64, dt: f64) -> Self {
        Self { t_i, t_f, dt }
    }

    /// Number of grid points, inclusive of both endpoints.
    pub fn len(&self) -> usize {
        ((self.t_f - self.t_i) / self.dt).round() as usize + 1
    }

    /// A grid always has at least one point; provided for lint-completeness.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The time [s] at grid index `idx`.
    pub fn time(&self, idx: usize) -> f64 {
        self.t_i + (idx as f64) * self.dt
    }

    /// Iterator over all grid times [s].
    pub fn times(&self) -> impl Iterator<Item = f64> + '_ {
        (0..self.len()).map(move |i| self.time(i))
    }
}

/// Tunable parameters for the three-step algorithm (Table III defaults).
#[derive(Debug, Clone)]
pub struct SolveParams {
    /// Coarse-sample count `|T^d|` for Algorithm 1 initialization (Table III: 20).
    pub n_coarse: usize,
    /// Initial candidate-time count `n_init` (Table III: 6).
    pub n_init: usize,
    /// Convergence tolerance `eps_cost` (Table III: 0.01).
    pub eps_cost: f64,
    /// Slack-removal tolerance `eps_remove` (Table III: 0.01).
    pub eps_remove: f64,
    /// Positive-definite weight `Q` for the extraction QP (identity in the example).
    pub q: SMatrix<f64, N, N>,
}

impl Default for SolveParams {
    fn default() -> Self {
        Self {
            n_coarse: 20,
            n_init: 6,
            eps_cost: 0.01,
            eps_remove: 0.01,
            q: SMatrix::<f64, N, N>::identity(),
        }
    }
}

/// The planner output: the maneuver set plus convergence diagnostics.
#[derive(Debug, Clone)]
pub struct Solution {
    /// Maneuvers `{t, Delta-v}`, one per optimal time in `T^opt`.
    pub maneuvers: Vec<Maneuver>,
    /// Total fuel cost sum of ||Delta-v_j|| [m/s].
    pub total_dv: f64,
    /// Algorithm 2 iteration count.
    pub iterations: usize,
    /// Relative residual ||w_err|| / ||w||.
    pub residual: f64,
    /// Optimal dual lambda_opt.
    pub lambda: Dual,
}

/// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1` for one candidate time.
///
/// Linear rows encode `a^T lambda <= b`; SOC rows encode `||G lambda||_2 <= h`.
/// [`crate::refine_socp`] assembles these into clarabel cones.
#[derive(Debug, Clone, Default)]
pub struct ConicRows {
    /// Linear rows `(a, b)` with `a^T lambda <= b`.
    pub linear: Vec<(SVector<f64, N>, f64)>,
    /// Second-order-cone rows `(G, h)` with `||G lambda||_2 <= h`.
    pub soc: Vec<(SMatrix<f64, M, N>, f64)>,
}

/// Primal fuel generator for one maneuver in the direct min-fuel SOCP
/// (Phase 5b, Algorithm 3). Describes how a Δv at one candidate time is built
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

/// Errors surfaced by the planner.
#[derive(Debug, Error)]
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
        /// Mean anomaly [rad].
        m: f64,
        /// Eccentricity.
        e: f64,
    },
    /// An input precondition was violated.
    #[error("invalid input: {0}")]
    InvalidInput(String),
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn dimensions_match_spec() {
        assert_eq!(N, 6);
        assert_eq!(M, 3);
    }

    #[test]
    fn worked_example_grid_has_3934_times() {
        let g = TimeGrid::uniform(0.0, 117990.0, 30.0);
        assert_eq!(g.len(), 3934);
        assert_abs_diff_eq!(g.time(0), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(g.time(g.len() - 1), 117990.0, epsilon = 1e-6);
        assert_eq!(g.times().count(), 3934);
    }

    #[test]
    fn hunter_grid_has_3901_times() {
        let g = TimeGrid::uniform(0.0, 39000.0, 10.0);
        assert_eq!(g.len(), 3901);
    }

    #[test]
    fn default_params_match_table_iii() {
        let p = SolveParams::default();
        assert_eq!(p.n_coarse, 20);
        assert_eq!(p.n_init, 6);
        assert_abs_diff_eq!(p.eps_cost, 0.01, epsilon = 1e-12);
        assert_abs_diff_eq!(p.eps_remove, 0.01, epsilon = 1e-12);
        assert_eq!(p.q, SMatrix::<f64, N, N>::identity());
    }

    #[test]
    fn conic_rows_default_is_empty() {
        let c = ConicRows::default();
        assert!(c.linear.is_empty());
        assert!(c.soc.is_empty());
    }
}
