//! Algorithm 3 - Control-Input Extraction: a direct min-fuel SOCP over the
//! converged active set `T^opt` recovers full 3-DOF maneuvers, including on
//! degenerate flat contacts where a fixed-direction recovery would under-span w.

use crate::cost::CostModel;
use crate::solver::min_fuel_socp;
use crate::types::{FuelGenerator, Maneuver, PlannerError, Pseudostate, TimeGrid, M, N};
use nalgebra::{SMatrix, SVector};

/// Maneuvers whose magnitude is below this fraction of the largest recovered
/// maneuver are interior-point dust and are pruned from the reported plan.
/// (Interior-point solutions are not exactly sparse, so a relative-magnitude
/// threshold is used to drop negligible maneuvers.)
const PRUNE_REL: f64 = 1e-3;

/// Relative tolerance for the primal/dual self-consistency gate: the extracted
/// min-fuel objective must agree with the refinement dual budget `c*` to within
/// this fraction. `budget` and this primal are the dual and primal of the **same**
/// converged active set `T^opt`, so conic strong duality ties them to interior-point
/// accuracy — ~1e-9 on the worked example (`total_dv` matches the refinement
/// objective to the printed digits), independent of `ε_cost`. `ε_cost` instead
/// bounds the *separate* gap between the restricted primal and the all-times dual
/// (0.62% here, the paper's 82.4-vs-82.0 mm/s sandwich), which this gate never sees.
/// So 5% is generous slack over interior-point noise on the ill-conditioned e=0.7
/// contacts while still catching a genuine primal/dual divergence.
const BUDGET_REL_TOL: f64 = 5e-2;

/// Result of Algorithm 3.
#[derive(Debug, Clone)]
pub(super) struct ExtractOutcome {
    pub maneuvers: Vec<Maneuver>,
    pub total_dv: f64,
    pub residual: f64,
}

/// Algorithm 3 — recover maneuvers over `T^opt` by direct min-fuel SOCP.
///
/// `budget` is the dual optimum `c*` from refinement (\[KD20\] eq. 40, Theorem 2).
/// It is not a constraint on the SOCP, but the extracted primal objective is
/// checked against it: by conic strong duality (\[KD20\] Theorems 1–3) the two
/// must agree, so a relative gap beyond [`BUDGET_REL_TOL`] is surfaced as
/// [`PlannerError::SolverFailed`]. This gate is **always on** (it was previously a
/// release-stripped `debug_assert!`).
///
/// Ref: \[KD20\] Algorithm 3 (Control Input Extraction); eq. 4 / eq. 7;
/// \[CD18\] eq. 2 / eq. 1.
///
/// # Errors
/// Propagates [`PlannerError::SolverFailed`] from the min-fuel SOCP, and returns
/// it when the recovered objective is inconsistent with `budget` (see above).
pub(super) fn extract<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    budget: f64,
    t_opt: &[usize],
) -> Result<ExtractOutcome, PlannerError> {
    // Build the per-candidate-time dynamics and fuel generators over T^opt.
    let gammas_t: Vec<SMatrix<f64, N, M>> = t_opt.iter().map(|&k| gammas[k]).collect();
    let generators: Vec<FuelGenerator> = t_opt
        .iter()
        .map(|&k| cost.at(grid.time(k)).fuel_generator())
        .collect();

    let sol = min_fuel_socp(w, &gammas_t, &generators)?;

    // Self-consistency gate (always on, not a release-stripped `debug_assert!`):
    // by conic strong duality the primal min-fuel objective must equal the
    // Algorithm-2 dual budget `c* = λ_optᵀw` (\[KD20\] Theorems 1–3; the eq.-37
    // lower bound). A gap beyond `BUDGET_REL_TOL` means the SOCP — which
    // `check_status` may accept at reduced accuracy (`AlmostSolved`) — returned a
    // primal inconsistent with the dual: a numerical failure the caller must see,
    // surfaced as `SolverFailed` rather than silently passing in release.
    //
    // The non-finite check is first and unconditional: a NaN/∞ objective must fail
    // the gate, but `NaN >= BUDGET_REL_TOL` is `false`, so the relative-gap test
    // alone would let it slip into `total_dv`. (`budget <= 0` only skips the
    // *relative* check — a degenerate non-positive dual carries no usable scale.)
    if !sol.objective.is_finite() {
        return Err(PlannerError::SolverFailed(format!(
            "min-fuel objective is non-finite ({})",
            sol.objective
        )));
    }
    if budget > 0.0 {
        let rel_gap = (sol.objective - budget).abs() / budget;
        if rel_gap >= BUDGET_REL_TOL {
            return Err(PlannerError::SolverFailed(format!(
                "min-fuel objective {} disagrees with dual budget {budget} \
                 (relative gap {rel_gap:.3} >= {BUDGET_REL_TOL})",
                sol.objective
            )));
        }
    }

    // Both `total_dv` and `residual` are reported on the FULL, pre-prune min-fuel
    // solution; only `maneuvers` is pruned of interior-point dust below.
    //
    // `total_dv` is the minimized fuel cost — the objective `Σⱼ f_{tⱼ}(Δvⱼ)` (the
    // paper's "delta-v cost" `c*`; eq. 4): `Σ‖Δvⱼ‖₂` under the L2 model, the
    // polytope gauge `Σθ` under FaceMax. It is the cost that was actually
    // minimized, NOT the L2 norm of the recovered net Δv (which under-states the
    // FaceMax gauge whenever a burn combines ≥2 vertices).
    let total_dv = sol.objective;

    // Residual of the FULL (unpruned) min-fuel solution: this is the ~0 we report.
    let mut w_acc_full = SVector::<f64, N>::zeros();
    for (idx, &k) in t_opt.iter().enumerate() {
        w_acc_full += gammas[k] * sol.dvs[idx];
    }
    let residual = (w - w_acc_full).norm() / w.norm();

    // Prune interior-point dust: keep maneuvers >= PRUNE_REL of the largest.
    let max_dv = sol.dvs.iter().map(|dv| dv.norm()).fold(0.0_f64, f64::max);
    let keep = PRUNE_REL * max_dv;
    let mut maneuvers = Vec::new();
    for (idx, &k) in t_opt.iter().enumerate() {
        let dv = sol.dvs[idx];
        if dv.norm() <= keep {
            continue;
        }
        maneuvers.push(Maneuver {
            t: grid.time(k),
            dv,
        });
    }

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
    use crate::types::PlannerError;

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
        fn gamma(&self, _t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
            Ok(top_identity())
        }
    }

    /// Γ(t) = top identity for t < 5, the zero matrix otherwise.
    struct TopThenZero;
    impl Dynamics for TopThenZero {
        fn gamma(&self, t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
            if t < 5.0 {
                Ok(top_identity())
            } else {
                Ok(SMatrix::<f64, N, M>::zeros())
            }
        }
    }

    #[test]
    fn extract_recovers_single_maneuver_with_zero_residual() {
        // w = (3,4,12,0,0,0); ‖w_top‖ = 13. min-fuel SOCP recovers dv = w_top
        // directly; budget = 13. Tolerance 1e-3: interior-point solvers achieve
        // ~5e-4 accuracy when the budget constraint is exactly tight.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0).unwrap();
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| TopId.gamma(t).unwrap()).collect();
        let cost = Piecewise::new(1.0e12).unwrap(); // Norm2
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let out = extract(&cost, &grid, &gammas, &w, 13.0, &[0]).unwrap();

        assert_eq!(out.maneuvers.len(), 1);
        let dv = out.maneuvers[0].dv;
        assert!((dv - SVector::<f64, M>::new(3.0, 4.0, 12.0)).norm() < 1e-3);
        assert!((out.total_dv - 13.0).abs() < 1e-3);
        assert!(out.residual < 1e-3);
    }

    #[test]
    fn extract_total_dv_is_the_facemax_gauge_not_the_l2_norm() {
        // Under FaceMax the reported `total_dv` must be
        // the minimized polytope gauge `Σθ` (the paper's "delta-v cost" `c*`,
        // eq. 4 / eq. 9), NOT the L2 norm of the recovered net Δv. The target
        // lands on a tetrahedron FACE, `w_top = v0 + v2`, whose gauge is 2.0
        // (θ0 = θ2 = 1) while `‖v0 + v2‖₂ = √(4/3) ≈ 1.1547`.
        let cols = crate::cost::facemax::vertex_columns();
        let v0 = cols[0];
        let v2 = cols[2];
        let target_dv = v0 + v2; // = (√(2/3), √(2/3), 0): ‖·‖₂ ≈ 1.1547, gauge = 2.0

        let grid = TimeGrid::uniform(0.0, 10.0, 1.0).unwrap();
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| TopId.gamma(t).unwrap()).collect();
        // Perigee at t = 0 ⇒ the t = 0 candidate is FaceMax-charged (Polytope gauge).
        let cost = Piecewise::with_perigee_epoch(40_000.0, 0.0).unwrap();
        let mut w = SVector::<f64, N>::zeros();
        for r in 0..M {
            w[r] = target_dv[r];
        }
        // budget = the true gauge optimum c* = 2.0.
        let out = extract(&cost, &grid, &gammas, &w, 2.0, &[0]).unwrap();

        assert_eq!(out.maneuvers.len(), 1);
        assert!((out.maneuvers[0].dv - target_dv).norm() < 1e-3);
        // The reported cost is the gauge (2.0), not the L2 net magnitude (1.1547).
        assert!(
            (out.total_dv - 2.0).abs() < 1e-3,
            "total_dv = {} (expected the FaceMax gauge 2.0, not the L2 norm 1.1547)",
            out.total_dv
        );
        assert!(out.residual < 1e-3);
    }

    #[test]
    fn extract_drops_unused_times() {
        // T^opt includes a time (t=8) where Γ = 0; min-fuel assigns it Δv≈0,
        // then pruning drops it. Only the t=0 maneuver survives.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0).unwrap();
        let gammas: Vec<SMatrix<f64, N, M>> = grid
            .times()
            .map(|t| TopThenZero.gamma(t).unwrap())
            .collect();
        let cost = Piecewise::new(1.0e12).unwrap();
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let out = extract(&cost, &grid, &gammas, &w, 13.0, &[0, 8]).unwrap();

        assert_eq!(out.maneuvers.len(), 1); // t=8 pruned
        assert!((out.maneuvers[0].t - 0.0).abs() < 1e-12);
        assert!(out.residual < 1e-3);
    }

    #[test]
    fn extract_rejects_objective_inconsistent_with_dual_budget() {
        // Self-consistency gate. By conic strong duality the extracted
        // min-fuel objective must equal the Algorithm-2 dual budget `c*`
        // ([KD20] Theorems 1-3). When it disagrees beyond tolerance, extraction
        // must surface `SolverFailed` — an always-on check, NOT the release-stripped
        // `debug_assert!` it used to be. Here the true min-fuel cost is 13, but a
        // wildly wrong budget of 100 is passed (relative gap 0.87 >> 5%).
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0).unwrap();
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| TopId.gamma(t).unwrap()).collect();
        let cost = Piecewise::new(1.0e12).unwrap(); // Norm2
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let err = extract(&cost, &grid, &gammas, &w, 100.0, &[0]).unwrap_err();
        assert!(
            matches!(err, PlannerError::SolverFailed(_)),
            "expected SolverFailed on a primal/dual mismatch, got {err:?}"
        );
    }

    #[test]
    fn extract_accepts_objective_matching_dual_budget() {
        // The gate must NOT fire on a consistent budget: true cost 13, budget 13.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0).unwrap();
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| TopId.gamma(t).unwrap()).collect();
        let cost = Piecewise::new(1.0e12).unwrap();
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        assert!(extract(&cost, &grid, &gammas, &w, 13.0, &[0]).is_ok());
    }
}
