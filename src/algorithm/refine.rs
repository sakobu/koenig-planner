//! Algorithm 2 - Iterative Refinement: solve eq. 40 on `T^est`, drop slack
//! times, add violated local maxima, until convergence.

use super::contact_on_grid;
use crate::cost::CostModel;
use crate::solver::refine_socp;
use crate::types::{ConicRows, Dual, PlannerError, Pseudostate, SolveParams, TimeGrid, M, N};
use nalgebra::SMatrix;

/// Values within `PLATEAU_EPS` of each other are treated as a flat top.
const PLATEAU_EPS: f64 = 1e-12;

/// Indices of `g` that are local maxima **and** exceed `threshold`.
///
/// A flat top (run of values within [`PLATEAU_EPS`]) yields a single
/// representative (the plateau midpoint). Endpoints are local maxima by the
/// boundary rule (compared only against their one in-bounds neighbor). The
/// global maximum is always a local maximum, so a violated global max is always
/// returned — guaranteeing Algorithm 2 makes progress.
pub(super) fn violated_local_maxima(g: &[f64], threshold: f64) -> Vec<usize> {
    let n = g.len();
    let mut out = Vec::new();
    let mut k = 0usize;
    while k < n {
        // Extent of the flat run [k..=j] of values ~equal to g[k].
        let mut j = k;
        while j + 1 < n && (g[j + 1] - g[k]).abs() <= PLATEAU_EPS {
            j += 1;
        }
        let left_ok = k == 0 || g[k - 1] < g[k];
        let right_ok = j == n - 1 || g[j + 1] < g[k];
        if left_ok && right_ok && g[k] > threshold {
            out.push((k + j) / 2); // plateau midpoint representative
        }
        k = j + 1;
    }
    out
}

/// Result of Algorithm 2.
#[derive(Debug, Clone)]
pub(super) struct RefineOutcome {
    /// Optimal candidate times `T^opt` (grid indices) — the active set that
    /// produced `lambda`.
    pub t_opt: Vec<usize>,
    /// Optimal dual `λ_opt`.
    pub lambda: Dual,
    /// Optimal objective `c* = λ_optᵀw` (the eq. 40 budget for extraction).
    pub objective: f64,
    /// Number of `refine_socp` solves performed.
    pub iterations: usize,
    /// `max_t g` after each solve — non-increasing; read only by tests.
    #[allow(dead_code)]
    pub max_g_trace: Vec<f64>,
}

/// Algorithm 2 — iteratively refine `T^est` until `max_t g ≤ 1 + ε_cost`.
pub(super) fn refine<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    params: &SolveParams,
    mut t_est: Vec<usize>,
    max_iters: usize,
) -> Result<RefineOutcome, PlannerError> {
    let add_threshold = 1.0;
    let converge_threshold = 1.0 + params.eps_cost;
    let keep_threshold = 1.0 - params.eps_remove;

    let mut max_g_trace = Vec::new();
    let mut iterations = 0usize;

    loop {
        if t_est.is_empty() {
            return Err(PlannerError::InvalidInput(
                "refine: candidate-time set became empty".into(),
            ));
        }

        // Solve eq. 40 over the current candidate set T^est.
        let rows: Vec<ConicRows> = t_est
            .iter()
            .map(|&k| cost.at(grid.time(k)).cone_constraints(&gammas[k]))
            .collect();
        let sol = refine_socp(w, &rows)?;
        iterations += 1;

        // Recompute g over the FULL grid with the new dual.
        let g = contact_on_grid(cost, grid, gammas, &sol.lambda);
        let max_g = g.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        max_g_trace.push(max_g);

        // Converged: T^opt is exactly the active set just solved.
        if max_g <= converge_threshold {
            let mut t_opt = t_est;
            t_opt.sort_unstable();
            t_opt.dedup();
            return Ok(RefineOutcome {
                t_opt,
                lambda: sol.lambda,
                objective: sol.objective,
                iterations,
                max_g_trace,
            });
        }
        if iterations >= max_iters {
            return Err(PlannerError::NotConverged {
                max_iters,
                achieved: max_g,
                target: converge_threshold,
            });
        }

        // Drop slack times, then add violated local maxima over the full grid.
        t_est.retain(|&k| g[k] >= keep_threshold);
        for k in violated_local_maxima(&g, add_threshold) {
            t_est.push(k);
        }
        t_est.sort_unstable();
        t_est.dedup();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::Piecewise;
    use crate::dynamics::Dynamics;
    use nalgebra::SVector;

    /// Mock dynamics: gently rotating, well-conditioned Γ(t). `rate` is chosen
    /// per test so directions vary across the grid without aliasing.
    struct SpinDyn {
        rate: f64,
    }
    impl Dynamics for SpinDyn {
        fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
            let a = self.rate * t;
            let (c, s) = (a.cos(), a.sin());
            SMatrix::<f64, N, M>::from_row_slice(&[
                c,
                -s,
                0.0, //
                s,
                c,
                0.0, //
                0.0,
                0.0,
                1.0, //
                0.5 * c,
                0.0,
                0.5 * s, //
                0.0,
                0.5 * c,
                -0.5 * s, //
                0.5 * s,
                -0.5 * c,
                0.0,
            ])
        }
    }

    fn cache(dynamics: &SpinDyn, grid: &TimeGrid) -> Vec<SMatrix<f64, N, M>> {
        grid.times().map(|t| dynamics.gamma(t)).collect()
    }

    #[test]
    fn refine_converges_and_trace_is_non_increasing() {
        let dynamics = SpinDyn { rate: 0.05 };
        let grid = TimeGrid::uniform(0.0, 60.0, 1.0); // 61 points
        let gammas = cache(&dynamics, &grid);
        let cost = Piecewise::new(1.0e12); // Norm2 everywhere
                                           // Reachable target: impulses at two distinct grid times.
        let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
        let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
        let w = dynamics.gamma(12.0) * ua + dynamics.gamma(47.0) * ub;
        let params = SolveParams::default();

        let out = refine(&cost, &grid, &gammas, &w, &params, vec![0, 30, 60], 50).unwrap();

        // Converged within tolerance.
        assert!(*out.max_g_trace.last().unwrap() <= 1.0 + params.eps_cost + 1e-9);
        assert!(!out.t_opt.is_empty());
        assert!(out.iterations >= 1 && out.iterations <= 50);
        // max_t g is monotonically non-increasing (small tol for solver noise).
        for pair in out.max_g_trace.windows(2) {
            assert!(
                pair[1] <= pair[0] + 1e-6,
                "trace not non-increasing: {:?}",
                out.max_g_trace
            );
        }
    }

    #[test]
    fn refine_reports_not_converged_at_iteration_cap() {
        let dynamics = SpinDyn { rate: 0.05 };
        let grid = TimeGrid::uniform(0.0, 60.0, 1.0);
        let gammas = cache(&dynamics, &grid);
        let cost = Piecewise::new(1.0e12);
        let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
        let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
        let w = dynamics.gamma(12.0) * ua + dynamics.gamma(47.0) * ub;
        let params = SolveParams::default();

        // Three candidate times don't span T^opt, so one solve won't converge.
        let err = refine(&cost, &grid, &gammas, &w, &params, vec![0, 30, 60], 1).unwrap_err();
        match err {
            PlannerError::NotConverged { max_iters, .. } => assert_eq!(max_iters, 1),
            other => panic!("expected NotConverged, got {other:?}"),
        }
    }

    #[test]
    fn refine_handles_mixed_facemax_and_norm2_cones() {
        // Risk R7: exercise the non-smooth FaceMax cost in the refine path.
        // Realistic period so the perigee window is a thin band, giving a true
        // FaceMax/Norm2 mix across the grid.
        let dynamics = SpinDyn { rate: 1.0e-4 };
        let grid = TimeGrid::uniform(0.0, 80_000.0, 1000.0); // 81 points
        let gammas = cache(&dynamics, &grid);
        let cost = Piecewise::new(40_000.0); // FaceMax for t mod 40000 in (16400, 23600)
                                             // Seed T^est with one Norm2 time (k=5, t=5000) and one FaceMax time
                                             // (k=20, t=20000) so the FIRST refine_socp assembles mixed cones.
        let ua = SVector::<f64, M>::new(0.5, -0.4, 0.6);
        // FaceMax support directions are tetrahedral vertices; use vertex 0
        // = [√(2/3), 0, -√(1/3)] so the FaceMax time is genuinely reachable.
        let v0 = SVector::<f64, M>::new((2.0_f64 / 3.0).sqrt(), 0.0, -(1.0_f64 / 3.0).sqrt());
        let w = dynamics.gamma(5000.0) * ua + dynamics.gamma(20_000.0) * (0.8 * v0);
        let params = SolveParams::default();

        // Mixed-cone refine_socp must succeed and the trace must be sane; we do
        // not hard-assert convergence-to-1 here (real conditioning is Phase 5's
        // job) — only that the non-smooth path runs and the dual is finite.
        let out = refine(&cost, &grid, &gammas, &w, &params, vec![5, 20], 50).unwrap();
        assert!(out.lambda.iter().all(|x| x.is_finite()));
        for pair in out.max_g_trace.windows(2) {
            assert!(pair[1] <= pair[0] + 1e-6);
        }
    }

    #[test]
    fn interior_peak() {
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn left_endpoint_peak() {
        assert_eq!(violated_local_maxima(&[3.0, 1.0, 0.5], 1.0), vec![0]);
    }

    #[test]
    fn right_endpoint_peak() {
        assert_eq!(violated_local_maxima(&[0.5, 1.0, 3.0], 1.0), vec![2]);
    }

    #[test]
    fn flat_top_two_yields_one_midpoint() {
        // plateau [1,2] -> midpoint (1+2)/2 = 1.
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 2.0, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn flat_top_three_yields_one_midpoint() {
        // plateau [1,3] -> midpoint (1+3)/2 = 2.
        assert_eq!(
            violated_local_maxima(&[0.0, 2.0, 2.0, 2.0, 0.0], 1.0),
            vec![2]
        );
    }

    #[test]
    fn monotone_increasing_picks_last() {
        assert_eq!(violated_local_maxima(&[0.0, 1.0, 2.0, 3.0], 1.0), vec![3]);
    }

    #[test]
    fn monotone_decreasing_picks_first() {
        assert_eq!(violated_local_maxima(&[3.0, 2.0, 1.0, 0.0], 1.0), vec![0]);
    }

    #[test]
    fn all_below_threshold_is_empty() {
        assert!(violated_local_maxima(&[0.1, 0.2, 0.1], 1.0).is_empty());
    }

    #[test]
    fn two_separated_peaks() {
        assert_eq!(
            violated_local_maxima(&[0.0, 2.0, 0.5, 3.0, 0.0], 1.0),
            vec![1, 3]
        );
    }

    #[test]
    fn threshold_filters_low_peak() {
        // peak at idx 1 (1.5 > 1) kept; peak at idx 3 (0.8 <= 1) dropped.
        assert_eq!(
            violated_local_maxima(&[0.0, 1.5, 0.0, 0.8, 0.0], 1.0),
            vec![1]
        );
    }

    #[test]
    fn single_element_above_and_below() {
        assert_eq!(violated_local_maxima(&[5.0], 1.0), vec![0]);
        assert!(violated_local_maxima(&[0.5], 1.0).is_empty());
    }
}
