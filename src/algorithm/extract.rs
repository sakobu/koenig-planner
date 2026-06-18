//! Algorithm 3 - Control-Input Extraction: optimal directions, then the QP for
//! magnitudes.

use crate::cost::CostModel;
use crate::solver::extract_qp;
use crate::types::{Dual, Maneuver, PlannerError, Pseudostate, TimeGrid, M, N};
use nalgebra::{SMatrix, SVector};

/// Below this support-direction norm a time contributes no usable maneuver
/// (a `y_j = 0` column leaves `α_j` unconstrained) and is dropped.
const SUPPORT_EPS: f64 = 1e-9;

/// Result of Algorithm 3.
#[derive(Debug, Clone)]
pub(super) struct ExtractOutcome {
    pub maneuvers: Vec<Maneuver>,
    pub total_dv: f64,
    pub residual: f64,
}

/// Algorithm 3 — recover maneuver directions `s_j` and magnitudes `α_j`.
#[allow(clippy::too_many_arguments)] // 8 params mirror Algorithm 3's inputs; collapsed into solve() in Task 5
pub(super) fn extract<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    q: &SMatrix<f64, N, N>,
    lambda: &Dual,
    budget: f64,
    t_opt: &[usize],
) -> Result<ExtractOutcome, PlannerError> {
    // Optimal directions, dropping ~zero-support times.
    let mut times: Vec<f64> = Vec::new();
    let mut dirs: Vec<SVector<f64, M>> = Vec::new();
    let mut ys: Vec<SVector<f64, N>> = Vec::new();
    for &k in t_opt {
        let y_dir = gammas[k].transpose() * lambda;
        let s = cost.at(grid.time(k)).support(y_dir);
        if s.norm() < SUPPORT_EPS {
            continue;
        }
        times.push(grid.time(k));
        ys.push(gammas[k] * s);
        dirs.push(s);
    }

    // QP for the nonnegative magnitudes (errors InvalidInput if ys is empty).
    let alpha = extract_qp(w, &ys, q, budget)?;

    let mut maneuvers = Vec::with_capacity(alpha.len());
    let mut w_acc = SVector::<f64, N>::zeros();
    let mut total_dv = 0.0;
    for (j, &a) in alpha.iter().enumerate() {
        let dv = a * dirs[j];
        total_dv += dv.norm();
        w_acc += a * ys[j];
        maneuvers.push(Maneuver { t: times[j], dv });
    }
    let residual = (w - w_acc).norm() / w.norm();

    Ok(ExtractOutcome {
        maneuvers,
        total_dv,
        residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::Piecewise;
    use crate::dynamics::Dynamics;

    fn top_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(i, i)] = 1.0;
        }
        g
    }

    /// Γ(t) = top 3×3 identity for every t.
    struct TopId;
    impl Dynamics for TopId {
        fn gamma(&self, _t: f64) -> SMatrix<f64, N, M> {
            top_identity()
        }
    }

    /// Γ(t) = top identity for t < 5, the zero matrix otherwise.
    struct TopThenZero;
    impl Dynamics for TopThenZero {
        fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
            if t < 5.0 {
                top_identity()
            } else {
                SMatrix::<f64, N, M>::zeros()
            }
        }
    }

    #[test]
    fn extract_recovers_single_maneuver_with_zero_residual() {
        // w = (3,4,12,0,0,0); ‖w_top‖ = 13. λ_opt ∥ w with unit top, so the
        // Norm2 support is unit(w_top); budget = λᵀw = 13; α = 13, dv = w_top.
        // Tolerance 1e-3: interior-point solvers achieve ~5e-4 accuracy when the
        // budget constraint is exactly tight (binding at the optimum).
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0);
        let gammas: Vec<SMatrix<f64, N, M>> = grid.times().map(|t| TopId.gamma(t)).collect();
        let cost = Piecewise::new(1.0e12); // Norm2
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let lambda = w / 13.0; // ∥ w, unit top half
        let q = SMatrix::<f64, N, N>::identity();
        let out = extract(&cost, &grid, &gammas, &w, &q, &lambda, 13.0, &[0]).unwrap();

        assert_eq!(out.maneuvers.len(), 1);
        let dv = out.maneuvers[0].dv;
        assert!((dv - SVector::<f64, M>::new(3.0, 4.0, 12.0)).norm() < 1e-3);
        assert!((out.total_dv - 13.0).abs() < 1e-3);
        assert!(out.residual < 1e-3);
    }

    #[test]
    fn extract_drops_zero_support_times() {
        // T^opt includes a time (t=8) where Γ = 0 -> support 0 -> dropped.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0);
        let gammas: Vec<SMatrix<f64, N, M>> = grid.times().map(|t| TopThenZero.gamma(t)).collect();
        let cost = Piecewise::new(1.0e12);
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let lambda = w / 13.0;
        let q = SMatrix::<f64, N, N>::identity();
        let out = extract(&cost, &grid, &gammas, &w, &q, &lambda, 13.0, &[0, 8]).unwrap();

        assert_eq!(out.maneuvers.len(), 1); // t=8 filtered out
        assert!((out.maneuvers[0].t - 0.0).abs() < 1e-12);
        assert!(out.residual < 1e-3);
    }
}
