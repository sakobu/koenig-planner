//! Field-exhaustive conversions between the core domain types and the api DTOs.
//!
//! Every conversion destructures its source with no `..`, so a new field in a
//! core type fails to compile here until it is handled — drift is impossible to
//! merge. `From` is used where it fits; error/param mappers are plain fns.

use crate::dto::{ApiError, SolveParamsDto};
use koenig_damico_planner::{PlannerError, SolveParams};

/// Map a [`PlannerError`] to a `"bad_request"` [`ApiError`].
pub(crate) fn bad_request(e: PlannerError) -> ApiError {
    ApiError {
        kind: "bad_request",
        message: e.to_string(),
    }
}

/// Map a [`PlannerError`] from the dispatch (solve/solve_from_initial_times)
/// to the correct [`ApiError`] kind.
///
/// `InvalidInput` is caller-fixable → `"bad_request"`.
/// All other variants indicate a numerically unsolvable / failed problem →
/// `"solver"`.  The match is exhaustive so a future new variant forces an
/// explicit decision here at compile time.
pub(crate) fn map_dispatch_error(e: koenig_damico_planner::PlannerError) -> ApiError {
    use koenig_damico_planner::PlannerError;
    match e {
        // Caller-fixable: bad request.
        PlannerError::InvalidInput(_) => ApiError {
            kind: "bad_request",
            message: e.to_string(),
        },
        // Well-formed request, numerically unsolvable / solver failure.
        PlannerError::SolverFailed(_)
        | PlannerError::NotConverged { .. }
        | PlannerError::KeplerDivergence { .. } => ApiError {
            kind: "solver",
            message: e.to_string(),
        },
    }
}

/// Merge optional overrides from [`SolveParamsDto`] into the Table III defaults.
pub(crate) fn resolve_params(dto: Option<SolveParamsDto>) -> SolveParams {
    let mut p = SolveParams::default();
    if let Some(d) = dto {
        if let Some(v) = d.n_coarse {
            p.n_coarse = v;
        }
        if let Some(v) = d.n_init {
            p.n_init = v;
        }
        if let Some(v) = d.eps_cost {
            p.eps_cost = v;
        }
        if let Some(v) = d.eps_remove {
            p.eps_remove = v;
        }
    }
    p
}
