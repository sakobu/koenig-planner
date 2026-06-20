//! Builds and solves eq. 40 over a candidate-time set (linear + SOC cones from
//! each time's `cone_constraints`), maximize -> minimize.

use crate::solver::{check_status, silent_settings};
use crate::types::{ConicRows, Dual, PlannerError, Pseudostate, M, N};
use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SupportedConeT,
};

/// The eq. 40 optimum over a candidate-time set: the dual `lambda` and the
/// optimal objective `c* = lambda^T w` (the minimum fuel cost).
///
/// The per-candidate-time contact value `g_{U(1,t)}(Gamma^T(t) lambda)` that
/// Algorithm 2 thresholds against is **not** returned here: it is recomputed by
/// the caller via [`crate::SublevelSet::contact`] (the caller scans `g` over the
/// full grid `T` anyway). Do not confuse it with clarabel's raw cone slack.
#[derive(Debug, Clone)]
pub struct RefineSolution {
    /// Optimal dual `lambda*` in R^6 (outward reachable-set normal).
    pub lambda: Dual,
    /// Optimal objective `c* = lambda*^T w` (>= 0).
    pub objective: f64,
}

/// Solve eq. 40 — `maximize lambda^T w s.t. g_{U(1,t)}(Gamma^T(t) lambda) <= 1`
/// for every candidate time — over `rows`, one [`ConicRows`] per candidate time.
///
/// `lambda` is a free (sign-unconstrained) variable: there is no cone on it.
/// Maps `maximize` to clarabel's `minimize` via `q = -w`; recovers
/// `c* = w . lambda` directly from the primal solution.
pub fn refine_socp(w: &Pseudostate, rows: &[ConicRows]) -> Result<RefineSolution, PlannerError> {
    let n_linear: usize = rows.iter().map(|r| r.linear.len()).sum();
    let n_soc: usize = rows.iter().map(|r| r.soc.len()).sum();
    let total_rows = n_linear + (M + 1) * n_soc;
    if total_rows == 0 {
        return Err(PlannerError::InvalidInput(
            "refine_socp: empty candidate-time set (objective is unbounded)".into(),
        ));
    }

    // Assemble A (total_rows x N) and b in lockstep with the cone vector:
    // all linear rows first (one NonnegativeCone), then one SOC block per SOC.
    let mut a_dense: Vec<Vec<f64>> = Vec::with_capacity(total_rows);
    let mut b: Vec<f64> = Vec::with_capacity(total_rows);
    let mut cones: Vec<SupportedConeT<f64>> = Vec::new();

    if n_linear > 0 {
        cones.push(NonnegativeConeT(n_linear));
        for cr in rows {
            for (a, bi) in &cr.linear {
                a_dense.push(a.iter().copied().collect());
                b.push(*bi);
            }
        }
    }
    for cr in rows {
        for (g, h) in &cr.soc {
            // Scalar-bound row: A all-zero, b = h  ->  s_0 = h.
            a_dense.push(vec![0.0; N]);
            b.push(*h);
            // Vector rows: A = -G, b = 0  ->  s_{1..M} = G lambda.
            for i in 0..M {
                a_dense.push((0..N).map(|c| -g[(i, c)]).collect());
                b.push(0.0);
            }
            cones.push(SecondOrderConeT(M + 1));
        }
    }

    let a_csc = CscMatrix::from(&a_dense);
    let p_csc = CscMatrix::<f64>::zeros((N, N));
    let q: Vec<f64> = (0..N).map(|i| -w[i]).collect();

    let mut solver = DefaultSolver::new(&p_csc, &q, &a_csc, &b, &cones, silent_settings())
        .map_err(|e| PlannerError::SolverFailed(format!("clarabel setup failed: {e:?}")))?;
    solver.solve();
    check_status(solver.solution.status)?;

    let lambda = Dual::from_iterator(solver.solution.x.iter().copied());
    let objective = w.dot(&lambda);
    Ok(RefineSolution { lambda, objective })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::{FaceMax, Norm2, SublevelSet};
    use crate::types::{ConicRows, M, N};
    use approx::assert_relative_eq;
    use nalgebra::{SMatrix, SVector};

    // Gamma (6x3) with top 3x3 = I, bottom 3x3 = 0. Then Gamma^T = [I_3 | 0],
    // so Gamma^T lambda = (lambda_1, lambda_2, lambda_3); Gamma v = [v; 0].
    fn gamma_top_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(i, i)] = 1.0;
        }
        g
    }

    // Gamma (6x3) with bottom 3x3 = I, top 3x3 = 0. Then Gamma^T = [0 | I_3],
    // so Gamma^T lambda = (lambda_4, lambda_5, lambda_6).
    fn gamma_bottom_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(M + i, i)] = 1.0;
        }
        g
    }

    fn w6(v: [f64; N]) -> SVector<f64, N> {
        SVector::<f64, N>::from_row_slice(&v)
    }

    // max_t g_{U(1,t)}(Gamma^T(t) lambda) recomputed straight from the rows,
    // normalized by each bound: dual feasibility requires this be <= 1 (+tol).
    fn max_contact(rows: &[ConicRows], lam: &SVector<f64, N>) -> f64 {
        let mut g = f64::NEG_INFINITY;
        for cr in rows {
            for (a, b) in &cr.linear {
                g = g.max(a.dot(lam) / b);
            }
            for (gmat, h) in &cr.soc {
                g = g.max((gmat * lam).norm() / h);
            }
        }
        g
    }

    #[test]
    fn empty_candidate_set_is_invalid_input() {
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let err = refine_socp(&w, &[]).unwrap_err();
        assert!(matches!(err, crate::types::PlannerError::InvalidInput(_)));
    }

    #[test]
    fn s1_single_soc_pins_three_coords() {
        // Norm2 at one time with Gamma^T = [I_3|0]; w hits only those 3 coords.
        // max (3,4,12).(l1,l2,l3) s.t. ||(l1,l2,l3)|| <= 1  ->  ||(3,4,12)|| = 13.
        let rows = vec![Norm2.cone_constraints(&gamma_top_identity())];
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        assert_relative_eq!(sol.objective, 13.0, epsilon = 1e-6);
        assert_relative_eq!(sol.objective, w.dot(&sol.lambda), epsilon = 1e-9);
        assert_relative_eq!(sol.lambda[0], 3.0 / 13.0, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[1], 4.0 / 13.0, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[2], 12.0 / 13.0, epsilon = 1e-6);
        // lambda_4..6 are a free direction (w is zero there): driven ~0.
        for i in 3..N {
            assert!(
                sol.lambda[i].abs() < 1e-4,
                "free coord {i} = {}",
                sol.lambda[i]
            );
        }
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn s2_face_max_lp_closed_form() {
        // FaceMax (LP path; no published closed-form reference) at one time,
        // Gamma=[I_3;0], w=(0,0,1,0,0,0). max l3 s.t. v_k.(l1,l2,l3) <= 1.
        // Binding v3,v4 give b*l3 <= 1 -> l3 = 1/b = sqrt(3); l2 pinned 0.
        let rows = vec![FaceMax.cone_constraints(&gamma_top_identity())];
        let w = w6([0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        let sqrt3 = 3.0_f64.sqrt();
        assert_relative_eq!(sol.objective, sqrt3, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[2], sqrt3, epsilon = 1e-6);
        assert!(sol.lambda[1].abs() < 1e-6);
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn s3_mixed_soc_and_lp_validates_cone_ordering() {
        // The realistic Piecewise case: one FaceMax time (4 linear rows on
        // l1..3) + one Norm2 time (SOC on l4..6). w=(0,0,1,0,0,1).
        // Separable: l3 = sqrt(3) (face-max) and l6 = 1 (||(l4,l5,l6)|| <= 1).
        let rows = vec![
            FaceMax.cone_constraints(&gamma_top_identity()),
            Norm2.cone_constraints(&gamma_bottom_identity()),
        ];
        let w = w6([0.0, 0.0, 1.0, 0.0, 0.0, 1.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        let sqrt3 = 3.0_f64.sqrt();
        assert_relative_eq!(sol.objective, sqrt3 + 1.0, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[2], sqrt3, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[5], 1.0, epsilon = 1e-6);
        assert!(sol.lambda[1].abs() < 1e-6);
        assert!(sol.lambda[3].abs() < 1e-4);
        assert!(sol.lambda[4].abs() < 1e-4);
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn objective_is_scale_equivariant() {
        // Scaling w by k>0 scales c* by k and leaves the lambda direction fixed.
        let rows = vec![Norm2.cone_constraints(&gamma_top_identity())];
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let base = refine_socp(&w, &rows).unwrap();
        let scaled = refine_socp(&(w * 2.5), &rows).unwrap();
        assert_relative_eq!(scaled.objective, 2.5 * base.objective, epsilon = 1e-6);
        assert_relative_eq!(scaled.lambda[0], base.lambda[0], epsilon = 1e-6);
    }

    #[test]
    fn zero_w_gives_zero_objective() {
        // w = 0 -> q = 0: a feasibility-only SOCP, optimum c* = w.lambda = 0.
        let rows = vec![Norm2.cone_constraints(&gamma_top_identity())];
        let w = w6([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        assert!(sol.objective.abs() < 1e-9);
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn unbounded_socp_maps_to_solver_failed() {
        // Norm2 on gamma_bottom_identity constrains only lambda_4..6; w points
        // along lambda_1 (unconstrained) -> objective unbounded -> clarabel
        // returns a non-Solved status, which must surface as SolverFailed
        // *through the wrapper* (proves the solve()->check_status wiring).
        let rows = vec![Norm2.cone_constraints(&gamma_bottom_identity())];
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let err = refine_socp(&w, &rows).unwrap_err();
        assert!(matches!(err, crate::types::PlannerError::SolverFailed(_)));
    }
}
