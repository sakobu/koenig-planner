//! Convex-solver wrappers around `clarabel`: the refinement SOCP (eq. 40,
//! `refine_socp`), the direct min-fuel SOCP that the live extraction path runs
//! (Algorithm 3, `min_fuel_socp`), and the legacy fixed-direction extraction QP
//! (`extract_qp`), plus shared settings/status helpers.

pub mod extract_qp;
pub mod min_fuel_socp;
pub mod refine_socp;

pub use extract_qp::extract_qp;
pub use min_fuel_socp::{min_fuel_socp, MinFuelSolution};
pub use refine_socp::{refine_socp, RefineSolution};

use crate::types::PlannerError;
use clarabel::solver::{DefaultSettings, DefaultSettingsBuilder, SolverStatus};

/// Default clarabel settings with logging suppressed (keeps the test/CI output
/// clean and the per-iteration SOCP solves quiet during refinement).
pub(crate) fn silent_settings() -> DefaultSettings<f64> {
    DefaultSettingsBuilder::default()
        .verbose(false)
        .build()
        .expect("default clarabel settings are always valid")
}

/// Map a clarabel terminal status to a planner result. `Solved` and
/// `AlmostSolved` (reduced accuracy) are accepted; every other status is a
/// failure whose message names the underlying clarabel status.
pub(crate) fn check_status(status: SolverStatus) -> Result<(), PlannerError> {
    match status {
        SolverStatus::Solved | SolverStatus::AlmostSolved => Ok(()),
        other => Err(PlannerError::SolverFailed(format!(
            "clarabel terminated with status {other:?}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clarabel::solver::SolverStatus;

    #[test]
    fn check_status_accepts_solved_and_almost_solved() {
        assert!(check_status(SolverStatus::Solved).is_ok());
        assert!(check_status(SolverStatus::AlmostSolved).is_ok());
    }

    #[test]
    fn check_status_rejects_failures_naming_the_status() {
        for bad in [
            SolverStatus::PrimalInfeasible,
            SolverStatus::DualInfeasible,
            SolverStatus::MaxIterations,
            SolverStatus::MaxTime,
            SolverStatus::NumericalError,
            SolverStatus::InsufficientProgress,
            SolverStatus::Unsolved,
        ] {
            let err = check_status(bad).unwrap_err();
            // The error message must name the underlying clarabel status,
            // so solver-failure debugging is not blind.
            assert!(format!("{err}").contains(&format!("{bad:?}")));
        }
    }

    #[test]
    fn silent_settings_are_non_verbose() {
        assert!(!silent_settings().verbose);
    }
}
