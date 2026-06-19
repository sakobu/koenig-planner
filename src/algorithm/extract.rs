//! Algorithm 3 - Control-Input Extraction: a direct min-fuel SOCP over the
//! converged active set `T^opt` recovers full 3-DOF maneuvers (robust on the
//! degenerate flat contacts where the fixed-support-direction QP under-spans w).

use crate::cost::CostModel;
use crate::solver::min_fuel_socp;
use crate::types::{FuelGenerator, Maneuver, PlannerError, Pseudostate, TimeGrid, M, N};
use nalgebra::{SMatrix, SVector};

/// Maneuvers whose magnitude is below this fraction of the largest recovered
/// maneuver are interior-point dust and are pruned from the reported plan.
/// (Interior-point solutions are not exactly sparse; see Decision D5.)
const PRUNE_REL: f64 = 1e-3;

/// Result of Algorithm 3.
#[derive(Debug, Clone)]
pub(super) struct ExtractOutcome {
    pub maneuvers: Vec<Maneuver>,
    pub total_dv: f64,
    pub residual: f64,
}

/// Algorithm 3 — recover maneuvers over `T^opt` by direct min-fuel SOCP.
///
/// `budget` is the dual optimum `c*` from refinement; it is used only as a
/// self-consistency sanity reference (the SOCP objective should match it to
/// solver tolerance), not as a constraint.
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
    debug_assert!(
        budget <= 0.0 || (sol.objective - budget).abs() / budget < 5e-2,
        "min-fuel objective {} disagrees with dual budget {budget}",
        sol.objective
    );

    // NOTE: residual is measured on the FULL unpruned solution (true reachability); total_dv below sums only the kept maneuvers. They are deliberately measured pre/post-prune.
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
    let mut total_dv = 0.0;
    for (idx, &k) in t_opt.iter().enumerate() {
        let dv = sol.dvs[idx];
        if dv.norm() <= keep {
            continue;
        }
        total_dv += dv.norm();
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
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0);
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| TopId.gamma(t).unwrap()).collect();
        let cost = Piecewise::new(1.0e12); // Norm2
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let out = extract(&cost, &grid, &gammas, &w, 13.0, &[0]).unwrap();

        assert_eq!(out.maneuvers.len(), 1);
        let dv = out.maneuvers[0].dv;
        assert!((dv - SVector::<f64, M>::new(3.0, 4.0, 12.0)).norm() < 1e-3);
        assert!((out.total_dv - 13.0).abs() < 1e-3);
        assert!(out.residual < 1e-3);
    }

    #[test]
    fn extract_drops_unused_times() {
        // T^opt includes a time (t=8) where Γ = 0; min-fuel assigns it Δv≈0,
        // then pruning drops it. Only the t=0 maneuver survives.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0);
        let gammas: Vec<SMatrix<f64, N, M>> = grid
            .times()
            .map(|t| TopThenZero.gamma(t).unwrap())
            .collect();
        let cost = Piecewise::new(1.0e12);
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let out = extract(&cost, &grid, &gammas, &w, 13.0, &[0, 8]).unwrap();

        assert_eq!(out.maneuvers.len(), 1); // t=8 pruned
        assert!((out.maneuvers[0].t - 0.0).abs() < 1e-12);
        assert!(out.residual < 1e-3);
    }
}
