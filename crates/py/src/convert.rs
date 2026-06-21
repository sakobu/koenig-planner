//! Conversions between the PyO3 surface and the `crates/api` DTOs.
//!
//! `From` cannot express the response conversion (building `Py<Maneuver>` needs
//! the GIL and is fallible), so it is a `Python`-taking fn. Both directions
//! destructure their api source with no `..`, so a new `crates/api` field breaks
//! compilation here until it is handled.

use crate::{Maneuver, Orbit, Solution};
use koenig_damico_planner_api::{
    ApiError, CostSpec, ManeuverDto, OrbitDto, SolveParamsDto, SolveRequest, SolveResponse,
};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// Map an [`ApiError`] to the matching Python exception.
pub(crate) fn api_err_to_py(e: ApiError) -> PyErr {
    match e.kind {
        "bad_request" => PyValueError::new_err(e.message),
        _ => PyRuntimeError::new_err(e.message),
    }
}

/// Assemble the api request from the Python-side inputs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_request(
    chief: &Orbit,
    t_i: f64,
    t_f: f64,
    dt: f64,
    w_metres: [f64; 6],
    cost: &str,
    period: Option<f64>,
    t_perigee0: Option<f64>,
    n_coarse: Option<usize>,
    n_init: Option<usize>,
    eps_cost: Option<f64>,
    eps_remove: Option<f64>,
    initial_times: Option<Vec<f64>>,
) -> PyResult<SolveRequest> {
    let cost = match cost {
        "norm2" => CostSpec::Norm2,
        "facemax" => CostSpec::FaceMax,
        "piecewise" => CostSpec::Piecewise { period, t_perigee0 },
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown cost {other:?}; expected one of \"norm2\", \"facemax\", \"piecewise\""
            )))
        }
    };
    let Orbit { a, e, i, raan, argp, mean_anom } = chief;
    Ok(SolveRequest {
        chief: OrbitDto {
            a: *a,
            e: *e,
            i: *i,
            raan: *raan,
            argp: *argp,
            mean_anom: *mean_anom,
        },
        t_i,
        t_f,
        dt,
        w_metres,
        cost,
        params: Some(SolveParamsDto { n_coarse, n_init, eps_cost, eps_remove }),
        initial_times,
    })
}

/// Convert an api [`SolveResponse`] into the Python `Solution` pyclass.
pub(crate) fn solution_to_py(py: Python<'_>, resp: SolveResponse) -> PyResult<Solution> {
    let SolveResponse { maneuvers, total_dv, iterations, residual, lambda } = resp;
    let maneuvers = maneuvers
        .iter()
        .map(|m| {
            let ManeuverDto { t, dv } = m;
            Py::new(py, Maneuver { t: *t, dv: (dv[0], dv[1], dv[2]) })
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(Solution { maneuvers, total_dv, iterations, residual, lambda: lambda.to_vec() })
}
