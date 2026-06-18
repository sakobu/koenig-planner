//! The Algorithm 3 QP: solve for nonnegative magnitudes that minimize the
//! weighted pseudostate residual.

use crate::solver::{check_status, silent_settings};
use crate::types::{PlannerError, Pseudostate, N};
use clarabel::algebra::CscMatrix;
use clarabel::solver::{DefaultSolver, IPSolver, NonnegativeConeT, SupportedConeT};
use nalgebra::{SMatrix, SVector};

/// Algorithm 3 QP: pick nonnegative magnitudes `alpha_j >= 0` (one per maneuver
/// direction) minimizing the `Q`-weighted residual `(w - sum_j alpha_j y_j)^T Q
/// (w - sum_j alpha_j y_j)`, subject to `sum_j alpha_j <= budget`
/// (`budget = lambda_opt^T w`).
///
/// `ys[j] = Gamma(t_j) . s_j` is the pseudostate contribution of the unit
/// support direction `s_j` at the j-th optimal time; the caller builds the
/// `Maneuver` as `dv = alpha_j . s_j` applied at `t_j`.
pub fn extract_qp(
    w: &Pseudostate,
    ys: &[SVector<f64, N>],
    q_weight: &SMatrix<f64, N, N>,
    budget: f64,
) -> Result<Vec<f64>, PlannerError> {
    let k = ys.len();
    if k == 0 {
        return Err(PlannerError::InvalidInput(
            "extract_qp: no maneuver directions".into(),
        ));
    }
    if budget < 0.0 {
        return Err(PlannerError::InvalidInput(format!(
            "extract_qp: budget must be non-negative, got {budget}"
        )));
    }

    // Symmetrize Q defensively so the triu(P) packing cannot drop an
    // asymmetric part. Q is PD (identity by default), so qsym is PD.
    let qsym = (q_weight + q_weight.transpose()) * 0.5;
    let qy: Vec<SVector<f64, N>> = ys.iter().map(|y| qsym * y).collect();
    let qw = qsym * w;

    // P = 2 Y^T Q Y, emitted UPPER-TRIANGULAR (strict-lower zeroed; CscMatrix
    // drops the zeros). Keep the full diagonal value (do not halve it).
    let mut p_dense: Vec<Vec<f64>> = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in i..k {
            p_dense[i][j] = 2.0 * ys[i].dot(&qy[j]);
        }
    }
    let p_csc = CscMatrix::from(&p_dense);

    // q = -2 Y^T Q w.
    let q: Vec<f64> = ys.iter().map(|y| -2.0 * y.dot(&qw)).collect();

    // A = [ -I_K ; 1^T ]  ((K+1) x K), b = [ 0_K ; budget ], one NonnegativeCone.
    let mut a_dense: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
    for (i, _) in ys.iter().enumerate() {
        let mut row = vec![0.0; k];
        row[i] = -1.0;
        a_dense.push(row);
    }
    a_dense.push(vec![1.0; k]);
    let a_csc = CscMatrix::from(&a_dense);

    let mut b = vec![0.0; k + 1];
    b[k] = budget;
    let cones: [SupportedConeT<f64>; 1] = [NonnegativeConeT(k + 1)];

    let mut solver = DefaultSolver::new(&p_csc, &q, &a_csc, &b, &cones, silent_settings())
        .map_err(|e| PlannerError::SolverFailed(format!("clarabel setup failed: {e:?}")))?;
    solver.solve();
    check_status(solver.solution.status)?;

    Ok(solver.solution.x.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::N;
    use approx::assert_relative_eq;
    use nalgebra::{SMatrix, SVector};

    fn e(i: usize) -> SVector<f64, N> {
        let mut v = SVector::<f64, N>::zeros();
        v[i] = 1.0;
        v
    }
    fn w6(v: [f64; N]) -> SVector<f64, N> {
        SVector::<f64, N>::from_row_slice(&v)
    }
    // Weighted residual (w - Y alpha)^T Q (w - Y alpha).
    fn weighted_obj(
        w: &SVector<f64, N>,
        ys: &[SVector<f64, N>],
        q: &SMatrix<f64, N, N>,
        alpha: &[f64],
    ) -> f64 {
        let mut werr = *w;
        for (a, y) in alpha.iter().zip(ys) {
            werr -= *a * *y;
        }
        (werr.transpose() * q * werr)[(0, 0)]
    }

    #[test]
    fn no_directions_is_invalid_input() {
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let q = SMatrix::<f64, N, N>::identity();
        let err = extract_qp(&w, &[], &q, 10.0).unwrap_err();
        assert!(matches!(err, crate::types::PlannerError::InvalidInput(_)));
    }

    #[test]
    fn qp_a_interior_optimum_exact_fit() {
        // y1=e1, w=2 e1, budget slack -> alpha=2, residual 0.
        let q = SMatrix::<f64, N, N>::identity();
        let a = extract_qp(&w6([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]), &[e(0)], &q, 10.0).unwrap();
        assert_eq!(a.len(), 1);
        assert_relative_eq!(a[0], 2.0, epsilon = 1e-6);
    }

    #[test]
    fn qp_b_budget_binds() {
        // Same as A but budget=1 -> alpha=1, residual 1.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let a = extract_qp(&w, &[e(0)], &q, 1.0).unwrap();
        assert_relative_eq!(a[0], 1.0, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &[e(0)], &q, &a), 1.0, epsilon = 1e-6);
    }

    #[test]
    fn qp_c_nonneg_binds() {
        // w = -3 e1, only nonneg direction available -> alpha=0, residual 3.
        let q = SMatrix::<f64, N, N>::identity();
        let a = extract_qp(&w6([-3.0, 0.0, 0.0, 0.0, 0.0, 0.0]), &[e(0)], &q, 10.0).unwrap();
        assert!(a[0].abs() < 1e-6);
    }

    #[test]
    fn qp_d_two_orthonormal_budget_binds() {
        // y1=e1,y2=e2, w=(2,3,..), budget=4 (sum unconstrained = 5).
        // Equal-shrink -> alpha=(1.5,2.5), weighted obj 0.5.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let ys = [e(0), e(1)];
        let a = extract_qp(&w, &ys, &q, 4.0).unwrap();
        assert_relative_eq!(a[0], 1.5, epsilon = 1e-6);
        assert_relative_eq!(a[1], 2.5, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &ys, &q, &a), 0.5, epsilon = 1e-6);
    }

    #[test]
    fn qp_e_weighted_q_budget_binds() {
        // Q = diag(1,4,1,1,1,1), y1=e1,y2=e2, w=(2,3,..), budget=4.
        // KKT: 2(2-a1)=8(3-a2) with a1+a2=4 -> alpha=(1.2,2.8), weighted obj 0.8.
        let mut q = SMatrix::<f64, N, N>::identity();
        q[(1, 1)] = 4.0;
        let w = w6([2.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let ys = [e(0), e(1)];
        let a = extract_qp(&w, &ys, &q, 4.0).unwrap();
        assert_relative_eq!(a[0], 1.2, epsilon = 1e-6);
        assert_relative_eq!(a[1], 2.8, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &ys, &q, &a), 0.8, epsilon = 1e-6);
    }

    #[test]
    fn qp_f_non_orthogonal_directions_exercise_off_diagonal_p() {
        // y1=e1, y2=e1+e2 (non-orthogonal -> P has off-diagonal terms).
        // w=(2,3,..): unconstrained LS gives alpha1<0, so nonneg binds ->
        // alpha=(0, 2.5), residual^2 = 0.5. Catches a triu-P packing bug.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let y2 = e(0) + e(1);
        let ys = [e(0), y2];
        let a = extract_qp(&w, &ys, &q, 10.0).unwrap();
        assert!(
            a[0].abs() < 1e-6,
            "alpha1 should hit the nonneg bound: {}",
            a[0]
        );
        assert_relative_eq!(a[1], 2.5, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &ys, &q, &a), 0.5, epsilon = 1e-6);
    }

    #[test]
    fn qp_residual_unique_when_p_singular() {
        // Duplicate directions -> Y^T Q Y singular -> alpha non-unique, but the
        // residual w - Y*alpha is unique. Assert on the residual, not alpha.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let ys = [e(0), e(0)];
        let a = extract_qp(&w, &ys, &q, 10.0).unwrap();
        assert_relative_eq!(a[0] + a[1], 2.0, epsilon = 1e-5);
        assert!(weighted_obj(&w, &ys, &q, &a) < 1e-8);
    }
}
