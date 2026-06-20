//! Algorithm 1 - Initialization: pick the `n_init` coarse times with the
//! largest contact value as `T^est`.

use super::contact_at;
use crate::cost::CostModel;
use crate::types::{Dual, SolveParams, TimeGrid, M, N};
use nalgebra::SMatrix;

/// Up to `n_coarse` evenly spaced grid indices (the coarse set `T^d`),
/// inclusive of both endpoints, clamped to `[1, grid_len]` and deduped.
///
/// Ref: \[KD20\] Algorithm 1 (builds the coarse candidate set `T^d`).
pub(super) fn coarse_indices(grid_len: usize, n_coarse: usize) -> Vec<usize> {
    debug_assert!(grid_len >= 1);
    let n = n_coarse.clamp(1, grid_len);
    if n == 1 {
        return vec![0];
    }
    let mut idx: Vec<usize> = (0..n)
        .map(|j| ((j as f64) * (grid_len - 1) as f64 / (n - 1) as f64).round() as usize)
        .collect();
    idx.dedup(); // idx is non-decreasing; drop any rounding collisions
    idx
}

/// The `n_init` coarse times with the largest contact `g_{U(1,t)}(Γᵀ(t)·lambda)`,
/// returned as sorted grid indices (`T^est`). `lambda` is the initial dual `∥ w`.
///
/// Ref: \[KD20\] Algorithm 1 (Initialization); eq. 30 / eq. 27 (contact g).
pub(super) fn initialize<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    lambda: &Dual,
    params: &SolveParams,
) -> Vec<usize> {
    let coarse = coarse_indices(grid.len(), params.n_coarse);
    let mut scored: Vec<(usize, f64)> = coarse
        .iter()
        .map(|&k| (k, contact_at(cost, grid, gammas, k, lambda)))
        .collect();
    // Largest contact first; tie-break by index for determinism.
    scored.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
    });
    let n_init = params.n_init.clamp(1, scored.len());
    let mut picked: Vec<usize> = scored.into_iter().take(n_init).map(|(k, _)| k).collect();
    picked.sort_unstable();
    picked
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::Piecewise;
    use crate::dynamics::Dynamics;
    use crate::types::PlannerError;
    use nalgebra::SVector;

    #[test]
    fn coarse_indices_span_endpoints_and_dedup() {
        assert_eq!(
            coarse_indices(101, 11),
            vec![0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        );
        assert_eq!(coarse_indices(1, 20), vec![0]); // clamped to grid_len
        assert_eq!(coarse_indices(5, 1), vec![0]); // n_coarse == 1
        let big = coarse_indices(3, 10); // clamp to 3 -> [0,1,2]
        assert_eq!(big, vec![0, 1, 2]);
    }

    /// Mock dynamics: Γ(t) is the top 3×3 identity scaled by (1 + t), so the
    /// Norm2 contact ‖Γᵀλ‖ grows monotonically with t.
    struct RampDyn;
    impl Dynamics for RampDyn {
        fn gamma(&self, t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
            let mut g = SMatrix::<f64, N, M>::zeros();
            let s = 1.0 + t;
            for i in 0..M {
                g[(i, i)] = s;
            }
            Ok(g)
        }
    }

    #[test]
    fn initialize_picks_largest_contact_times() {
        let grid = TimeGrid::uniform(0.0, 100.0, 1.0).unwrap(); // 101 points
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| RampDyn.gamma(t).unwrap()).collect();
        let cost = Piecewise::new(1.0e12); // huge period -> Norm2 everywhere
        let w = SVector::<f64, N>::from_row_slice(&[1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);
        let params = SolveParams {
            n_coarse: 11,
            n_init: 3,
            ..SolveParams::default()
        };
        // Coarse set = {0,10,...,100}; contact grows with t, so the 3 largest
        // are the last three coarse indices, returned sorted.
        let t_est = initialize(&cost, &grid, &gammas, &w, &params);
        assert_eq!(t_est, vec![80, 90, 100]);
    }
}
