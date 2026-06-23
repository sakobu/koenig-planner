//! Field-exhaustive conversions between the core domain types and the api DTOs.
//!
//! Every conversion destructures its source with no `..`, so a new field in a
//! core type fails to compile here until it is handled — drift is impossible to
//! merge. `From` is used where it fits; error/param mappers are plain fns.

use crate::dto::{ApiError, ApiErrorKind, ManeuverDto, SolveParamsDto, SolveResponse};
use koenig_damico_planner::{Maneuver, PlannerError, PrimerHistory, Solution, SolveParams};

/// Map a [`PlannerError`] to a `"bad_request"` [`ApiError`].
pub(crate) fn bad_request(e: PlannerError) -> ApiError {
    ApiError {
        kind: ApiErrorKind::BadRequest,
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
            kind: ApiErrorKind::BadRequest,
            message: e.to_string(),
        },
        // Well-formed request, numerically unsolvable / solver failure.
        PlannerError::SolverFailed(_)
        | PlannerError::NotConverged { .. }
        | PlannerError::KeplerDivergence { .. } => ApiError {
            kind: ApiErrorKind::Solver,
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

impl From<&Maneuver> for ManeuverDto {
    fn from(m: &Maneuver) -> Self {
        let Maneuver { t, dv } = *m; // Maneuver is Copy
        ManeuverDto {
            t,
            dv: [dv[0], dv[1], dv[2]],
        }
    }
}

/// Build the wire response from the core [`Solution`] plus its
/// [`PrimerHistory`] (computed in `run` from the converged dual). Both sources
/// are destructured with no `..`, so a new field on either forces a decision
/// here at compile time.
impl From<(Solution, PrimerHistory)> for SolveResponse {
    fn from((sol, primer): (Solution, PrimerHistory)) -> Self {
        let Solution {
            maneuvers,
            total_dv,
            iterations,
            residual,
            lambda,
        } = sol;
        let PrimerHistory {
            times,
            vectors,
            magnitudes,
        } = primer;
        SolveResponse {
            maneuvers: maneuvers.iter().map(ManeuverDto::from).collect(),
            total_dv,
            iterations,
            residual,
            lambda: [
                lambda[0], lambda[1], lambda[2], lambda[3], lambda[4], lambda[5],
            ],
            primer_times: times,
            primer_magnitude: magnitudes,
            primer_rtn: vectors.iter().map(|p| [p[0], p[1], p[2]]).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::dto::SolveResponse;
    use koenig_damico_planner::{Maneuver, PrimerHistory, Solution};
    use nalgebra::{SVector, Vector3};

    #[test]
    fn solution_converts_field_for_field() {
        let sol = Solution {
            maneuvers: vec![Maneuver {
                t: 12.0,
                dv: Vector3::new(1.0, 2.0, 3.0),
            }],
            total_dv: 6.0,
            iterations: 4,
            residual: 1e-12,
            lambda: SVector::<f64, 6>::from_row_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        };
        let primer = PrimerHistory {
            times: vec![0.0, 30.0],
            vectors: vec![Vector3::new(0.1, 0.2, 0.3), Vector3::new(0.4, 0.5, 0.6)],
            magnitudes: vec![0.5, 1.0],
        };
        let resp: SolveResponse = (sol, primer).into();
        assert_eq!(resp.maneuvers.len(), 1);
        assert_eq!(resp.maneuvers[0].t, 12.0);
        assert_eq!(resp.maneuvers[0].dv, [1.0, 2.0, 3.0]);
        assert_eq!(resp.total_dv, 6.0);
        assert_eq!(resp.iterations, 4);
        assert_eq!(resp.residual, 1e-12);
        assert_eq!(resp.lambda, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(resp.primer_times, [0.0, 30.0]);
        assert_eq!(resp.primer_magnitude, [0.5, 1.0]);
        assert_eq!(resp.primer_rtn, [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]);
    }
}
