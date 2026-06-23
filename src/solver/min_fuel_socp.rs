//! Direct min-fuel SOCP: recover full 3-DOF maneuvers over a fixed
//! candidate-time set by minimizing the (piecewise) fuel cost subject to exact
//! reachability `Σⱼ Γ(tⱼ)·Δvⱼ = w`.
//!
//! Unlike the fixed-support-direction magnitude QP (`extract_qp`), this is robust
//! on the degenerate flat contacts of e=0.7 orbits, where the per-time support
//! directions are nearly collinear and cannot span `w`. Freeing each maneuver to
//! a full vector (charged by its true cost) recovers `w` to ~0 residual; by conic
//! strong duality the optimum equals the eq.40 dual value `c*`. Sum-of-norms is
//! group-sparse, so few maneuvers come out nonzero.

use crate::solver::{check_status, silent_settings};
use crate::types::{FuelGenerator, InvalidInputKind, PlannerError, Pseudostate, M, N};
use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SupportedConeT, ZeroConeT,
};
use nalgebra::{SMatrix, SVector};

/// Result of the direct min-fuel SOCP.
///
/// Ref: \[KD20\] eq. 42 (the recovered per-time maneuvers reaching `w`).
#[derive(Debug, Clone)]
pub struct MinFuelSolution {
    /// Recovered Δv per input candidate time (aligned to `gammas`/`generators`;
    /// times the optimum does not use come back ≈ 0 via the sum-of-norms penalty).
    pub dvs: Vec<SVector<f64, M>>,
    /// Total fuel cost `Σⱼ f_{tⱼ}(Δvⱼ)` (≈ the dual budget `c*`).
    pub objective: f64,
}

/// Solve `min Σⱼ f_{tⱼ}(Δvⱼ) s.t. Σⱼ Γ(tⱼ)·Δvⱼ = w` over the candidate times
/// described pairwise by (`gammas`, `generators`).
///
/// `f_{tⱼ}` is the cost whose unit ball the generator describes: an L2 norm
/// (`Norm` → one second-order cone) or a polytopic gauge (`Polytope` →
/// nonnegative-cone LP over the vertex directions). Returns one Δv per candidate
/// time plus the optimal fuel.
///
/// Ref: \[KD20\] eq. 4 (master min-fuel objective); eq. 33; eq. 9; \[CD18\] eq. 2
/// (sum-of-norms).
///
/// # Errors
/// - [`PlannerError::InvalidInput`] if `gammas` is empty or its length does not
///   match `generators` (one generator per candidate time is required).
/// - [`PlannerError::SolverFailed`] if the clarabel SOCP fails to set up or does
///   not reach a (near-)optimal status.
pub fn min_fuel_socp(
    w: &Pseudostate,
    gammas: &[SMatrix<f64, N, M>],
    generators: &[FuelGenerator],
) -> Result<MinFuelSolution, PlannerError> {
    let k = gammas.len();
    if k == 0 || k != generators.len() {
        return Err(PlannerError::InvalidInput(
            InvalidInputKind::EmptyCandidateSet,
        ));
    }

    // --- Per-maneuver variable layout (contiguous block per maneuver). ---
    //   Norm        -> [c_j, v_j(M)]   (M+1 vars; Δv_j = v_j; cost = c_j ≥ ‖v_j‖)
    //   Polytope(p) -> [θ_j(p)]        (p   vars; Δv_j = Σₖ θ_jk d_k; cost = Σ θ)
    let mut offsets = Vec::with_capacity(k);
    let mut n_var = 0usize;
    for generator in generators {
        offsets.push(n_var);
        n_var += match generator {
            FuelGenerator::Norm => M + 1,
            FuelGenerator::Polytope(dirs) => dirs.len(),
        };
    }

    // --- Objective q (P = 0): minimize the sum of the cost variables. ---
    let mut q = vec![0.0f64; n_var];
    for (j, generator) in generators.iter().enumerate() {
        match generator {
            FuelGenerator::Norm => q[offsets[j]] = 1.0, // the epigraph c_j
            FuelGenerator::Polytope(dirs) => {
                for kk in 0..dirs.len() {
                    q[offsets[j] + kk] = 1.0;
                }
            }
        }
    }

    // --- Equality-constraint column vectors per variable: the 6-vector that
    //     multiplies each variable in `Σ Γ Δv = w`. For Norm, the c_j column is
    //     0 and the v_j columns are Γ[:,m]; for Polytope, column k is Γ·d_k. ---
    let mut eq_cols: Vec<Vec<SVector<f64, N>>> = Vec::with_capacity(k);
    for (j, generator) in generators.iter().enumerate() {
        let g = &gammas[j];
        match generator {
            FuelGenerator::Norm => {
                let mut cols = vec![SVector::<f64, N>::zeros()]; // c_j
                for m in 0..M {
                    cols.push(SVector::<f64, N>::from_iterator(
                        g.column(m).iter().copied(),
                    ));
                }
                eq_cols.push(cols);
            }
            FuelGenerator::Polytope(dirs) => {
                eq_cols.push(dirs.iter().map(|d| g * d).collect());
            }
        }
    }

    // --- Assemble A, b, cones in lockstep:
    //     (1) N equality rows  -> ZeroCone(N)
    //     (2) per maneuver:  Norm -> SOC(M+1) on (c_j, v_j);  Polytope -> Nonneg(p)
    let mut a_rows: Vec<Vec<f64>> = Vec::new();
    let mut b: Vec<f64> = Vec::new();
    let mut cones: Vec<SupportedConeT<f64>> = Vec::new();

    // (1) Equality block: A·x = w  (s ∈ ZeroCone).
    for r in 0..N {
        let mut row = vec![0.0; n_var];
        for j in 0..k {
            for (c, col) in eq_cols[j].iter().enumerate() {
                row[offsets[j] + c] = col[r];
            }
        }
        a_rows.push(row);
        b.push(w[r]);
    }
    cones.push(ZeroConeT(N));

    // (2) Per-maneuver cost cones.  clarabel: s = b - A·x.
    for (j, generator) in generators.iter().enumerate() {
        match generator {
            FuelGenerator::Norm => {
                // s = (c_j, v_j) ∈ SOC(M+1): A = -selector, b = 0.
                let mut s0 = vec![0.0; n_var];
                s0[offsets[j]] = -1.0; // s_0 = c_j
                a_rows.push(s0);
                b.push(0.0);
                for m in 0..M {
                    let mut row = vec![0.0; n_var];
                    row[offsets[j] + 1 + m] = -1.0; // s_{1+m} = v_{j,m}
                    a_rows.push(row);
                    b.push(0.0);
                }
                cones.push(SecondOrderConeT(M + 1));
            }
            FuelGenerator::Polytope(dirs) => {
                // s_k = θ_jk ≥ 0.
                for kk in 0..dirs.len() {
                    let mut row = vec![0.0; n_var];
                    row[offsets[j] + kk] = -1.0;
                    a_rows.push(row);
                    b.push(0.0);
                }
                cones.push(NonnegativeConeT(dirs.len()));
            }
        }
    }

    let a_csc = CscMatrix::from(&a_rows);
    let p_csc = CscMatrix::<f64>::zeros((n_var, n_var));

    let mut solver = DefaultSolver::new(&p_csc, &q, &a_csc, &b, &cones, silent_settings())
        .map_err(|e| PlannerError::SolverFailed(format!("clarabel setup failed: {e:?}")))?;
    solver.solve();
    check_status(solver.solution.status)?;
    let x = &solver.solution.x;

    // --- Reconstruct Δv per maneuver and total fuel. ---
    let mut dvs = Vec::with_capacity(k);
    let mut objective = 0.0;
    for (j, generator) in generators.iter().enumerate() {
        match generator {
            FuelGenerator::Norm => {
                let v = SVector::<f64, M>::from_iterator((0..M).map(|m| x[offsets[j] + 1 + m]));
                objective += x[offsets[j]]; // c_j ≈ ‖v‖
                dvs.push(v);
            }
            FuelGenerator::Polytope(dirs) => {
                let mut v = SVector::<f64, M>::zeros();
                for (kk, d) in dirs.iter().enumerate() {
                    let theta = x[offsets[j] + kk];
                    v += theta * d;
                    objective += theta;
                }
                dvs.push(v);
            }
        }
    }

    Ok(MinFuelSolution { dvs, objective })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Γ (6x3) with top 3x3 = I, bottom 3x3 = 0: Γv = (v; 0), Γᵀλ = (λ₁,λ₂,λ₃).
    fn gamma_top_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(i, i)] = 1.0;
        }
        g
    }

    fn norm_gen() -> FuelGenerator {
        FuelGenerator::Norm
    }

    // The four FaceMax V_vertex columns, shared with the cost module so the two
    // cannot drift. Ref: [KD20] eq. 47.
    fn vertex_dirs() -> Vec<SVector<f64, M>> {
        crate::cost::facemax::vertex_columns().to_vec()
    }

    fn w6(v: [f64; N]) -> Pseudostate {
        Pseudostate::from_row_slice(&v)
    }

    // Ref: [KD20] eq. 4.
    fn residual(w: &Pseudostate, gammas: &[SMatrix<f64, N, M>], dvs: &[SVector<f64, M>]) -> f64 {
        let mut acc = SVector::<f64, N>::zeros();
        for (g, dv) in gammas.iter().zip(dvs) {
            acc += g * dv;
        }
        (w - acc).norm() / w.norm()
    }

    #[test]
    fn empty_or_mismatched_input_is_invalid() {
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(matches!(
            min_fuel_socp(&w, &[], &[]).unwrap_err(),
            PlannerError::InvalidInput(InvalidInputKind::EmptyCandidateSet)
        ));
        assert!(matches!(
            min_fuel_socp(&w, &[gamma_top_identity()], &[]).unwrap_err(),
            PlannerError::InvalidInput(InvalidInputKind::EmptyCandidateSet)
        ));
    }

    // Ref: [KD20] eq. 4.
    #[test]
    fn single_norm_time_recovers_exact_dv() {
        // w = (3,4,12,0,0,0), Γ = [I;0]: Δv = (3,4,12), fuel = ‖Δv‖ = 13.
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let gammas = [gamma_top_identity()];
        let sol = min_fuel_socp(&w, &gammas, &[norm_gen()]).unwrap();
        assert_relative_eq!(
            sol.dvs[0],
            SVector::<f64, M>::new(3.0, 4.0, 12.0),
            epsilon = 1e-5
        );
        assert_relative_eq!(sol.objective, 13.0, epsilon = 1e-4);
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
    }

    // Ref: [KD20] eq. 4.
    #[test]
    fn degenerate_collinear_support_still_spans_w() {
        // The degenerate-collinear scenario in miniature. Two times whose Γ map
        // control space to DIFFERENT pseudostate planes; the per-time support directions of a
        // shared λ are collinear (the failure mode of extract_qp), yet w needs a
        // contribution from each time. Full-DOF min-fuel recovers w exactly.
        // Time A: Γ_A = [I;0] (controls coords 0..2). Time B: Γ_B maps to coords 3..5.
        let mut gb = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            gb[(M + i, i)] = 1.0; // bottom identity
        }
        let gammas = [gamma_top_identity(), gb];
        let w = w6([1.0, 0.0, 0.0, 0.0, 2.0, 0.0]); // needs BOTH times
        let sol = min_fuel_socp(&w, &gammas, &[norm_gen(), norm_gen()]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5, "residual too large");
        // Fuel = ‖(1,0,0)‖ + ‖(0,2,0)‖ = 1 + 2 = 3.
        assert_relative_eq!(sol.objective, 3.0, epsilon = 1e-4);
    }

    #[test]
    fn unused_time_comes_back_zero() {
        // Second time contributes nothing to a reachable w -> its Δv ≈ 0.
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let gammas = [gamma_top_identity(), gamma_top_identity()];
        let sol = min_fuel_socp(&w, &gammas, &[norm_gen(), norm_gen()]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        // Total fuel is still 13 (split across the two identical times any way).
        assert_relative_eq!(sol.objective, 13.0, epsilon = 1e-3);
    }

    // Ref: [KD20] eq. 9; eq. 47.
    #[test]
    fn facemax_single_vertex_costs_unit() {
        // w = Γ·v0 with Γ = [I;0], v0 the first vertex direction (a unit vector).
        // The gauge of a vertex is 1, so fuel = 1 and Δv ≈ v0.
        let dirs = vertex_dirs();
        let gammas = [gamma_top_identity()];
        let mut w = SVector::<f64, N>::zeros();
        for r in 0..M {
            w[r] = dirs[0][r];
        }
        let sol = min_fuel_socp(&w, &gammas, &[FuelGenerator::Polytope(dirs.clone())]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        assert_relative_eq!(sol.objective, 1.0, epsilon = 1e-3);
        assert_relative_eq!(sol.dvs[0], dirs[0], epsilon = 1e-3);
    }

    // Ref: [KD20] eq. 49.
    #[test]
    fn mixed_norm_and_facemax_times() {
        // One Norm time (coords 0..2) + one FaceMax time (coords 3..5). Both used.
        let dirs = vertex_dirs();
        let mut gb = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            gb[(M + i, i)] = 1.0;
        }
        let gammas = [gamma_top_identity(), gb];
        // w: (1,0,0) on the Norm side; gb·v2 (a vertex) on the FaceMax side.
        let mut w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let face_part = gb * dirs[2];
        for r in 0..N {
            w[r] += face_part[r];
        }
        let sol = min_fuel_socp(
            &w,
            &gammas,
            &[norm_gen(), FuelGenerator::Polytope(dirs.clone())],
        )
        .unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        // Fuel = ‖(1,0,0)‖ (Norm) + 1 (vertex gauge) = 2.
        assert_relative_eq!(sol.objective, 2.0, epsilon = 1e-3);
    }

    // Ref: [KD20] eq. 9; eq. 48.
    #[test]
    fn facemax_two_vertex_combination_charges_sum() {
        // Target on a FACE, not a vertex: w = Gamma . (v0 + v2). The gauge must
        // COMBINE two vertices -> f(v0+v2) = 2 (theta0 = theta2 = 1; no cheaper
        // nonnegative combo exists, since the four vertices sum to 0). Exercises
        // the Polytope LP's multi-vertex path (single-vertex tests cannot).
        let dirs = vertex_dirs();
        let gammas = [gamma_top_identity()];
        let target_dv = dirs[0] + dirs[2];
        let mut w = SVector::<f64, N>::zeros();
        for r in 0..M {
            w[r] = target_dv[r];
        }
        let sol = min_fuel_socp(&w, &gammas, &[FuelGenerator::Polytope(dirs.clone())]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        assert_relative_eq!(sol.objective, 2.0, epsilon = 1e-3);
        assert_relative_eq!(sol.dvs[0], target_dv, epsilon = 1e-3);
    }

    #[test]
    fn unreachable_w_is_solver_failed() {
        // Γ = [I;0] reaches only coords 0..2; w in coord 4 is unreachable ->
        // the equality is primal-infeasible -> SolverFailed through the wrapper.
        let w = w6([0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
        let gammas = [gamma_top_identity()];
        let err = min_fuel_socp(&w, &gammas, &[norm_gen()]).unwrap_err();
        assert!(matches!(err, PlannerError::SolverFailed(_)));
    }
}
