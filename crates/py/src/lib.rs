//! # koenig_planner (Python bindings)
//!
//! Thin PyO3 wrapper over [`koenig_planner_api`]. See `crates/py/tests` for the
//! golden worked-example checks and `crates/py/examples` for a plotting showcase.

use koenig_planner_api::{run, ApiError, CostSpec, OrbitDto, SolveParamsDto, SolveRequest};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// Map an [`ApiError`] to the matching Python exception.
fn api_err_to_py(e: ApiError) -> PyErr {
    match e.kind {
        "bad_request" => PyValueError::new_err(e.message),
        _ => PyRuntimeError::new_err(e.message),
    }
}

/// Chief mean absolute orbit. Angles in **degrees**; `a` in **metres**.
#[pyclass(from_py_object)]
#[derive(Clone)]
struct Orbit {
    #[pyo3(get)]
    a: f64,
    #[pyo3(get)]
    e: f64,
    #[pyo3(get)]
    i: f64,
    #[pyo3(get)]
    raan: f64,
    #[pyo3(get)]
    argp: f64,
    #[pyo3(get)]
    mean_anom: f64,
}

#[pymethods]
impl Orbit {
    #[new]
    fn new(a: f64, e: f64, i: f64, raan: f64, argp: f64, mean_anom: f64) -> Self {
        Self {
            a,
            e,
            i,
            raan,
            argp,
            mean_anom,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Orbit(a={}, e={}, i={}, raan={}, argp={}, mean_anom={})",
            self.a, self.e, self.i, self.raan, self.argp, self.mean_anom
        )
    }
}

/// One impulsive maneuver: time `t` `[s]` and RTN delta-v `[m/s]`.
#[pyclass]
struct Maneuver {
    #[pyo3(get)]
    t: f64,
    #[pyo3(get)]
    dv: (f64, f64, f64),
}

#[pymethods]
impl Maneuver {
    fn __repr__(&self) -> String {
        format!(
            "Maneuver(t={}, dv=({}, {}, {}))",
            self.t, self.dv.0, self.dv.1, self.dv.2
        )
    }
}

/// Planner output.
#[pyclass]
struct Solution {
    #[pyo3(get)]
    maneuvers: Vec<Py<Maneuver>>,
    #[pyo3(get)]
    total_dv: f64,
    #[pyo3(get)]
    iterations: usize,
    #[pyo3(get)]
    residual: f64,
    /// Optimal dual `λ_opt ∈ ℝ⁶`. Named `lambda_` (trailing underscore) because
    /// `lambda` is a Python keyword.
    #[pyo3(get, name = "lambda_")]
    lambda: Vec<f64>,
}

#[pymethods]
impl Solution {
    fn __repr__(&self) -> String {
        format!(
            "Solution(maneuvers={}, total_dv={}, iterations={}, residual={})",
            self.maneuvers.len(),
            self.total_dv,
            self.iterations,
            self.residual
        )
    }
}

/// Plan a maneuver set.
///
/// `cost` is one of `"norm2"`, `"facemax"`, `"piecewise"`. `period` /
/// `t_perigee0` apply only to `"piecewise"` (defaults derived from the chief).
#[pyfunction]
#[pyo3(signature = (
    chief, t_i, t_f, dt, w_metres, cost="piecewise",
    *, period=None, t_perigee0=None,
    n_coarse=None, n_init=None, eps_cost=None, eps_remove=None,
    initial_times=None
))]
#[allow(clippy::too_many_arguments)]
fn solve(
    py: Python<'_>,
    chief: Orbit,
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
) -> PyResult<Solution> {
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

    let req = SolveRequest {
        chief: OrbitDto {
            a: chief.a,
            e: chief.e,
            i: chief.i,
            raan: chief.raan,
            argp: chief.argp,
            mean_anom: chief.mean_anom,
        },
        t_i,
        t_f,
        dt,
        w_metres,
        cost,
        params: Some(SolveParamsDto {
            n_coarse,
            n_init,
            eps_cost,
            eps_remove,
        }),
        initial_times,
    };

    let resp = run(req).map_err(api_err_to_py)?;

    let maneuvers = resp
        .maneuvers
        .iter()
        .map(|m| {
            Py::new(
                py,
                Maneuver {
                    t: m.t,
                    dv: (m.dv[0], m.dv[1], m.dv[2]),
                },
            )
        })
        .collect::<PyResult<Vec<_>>>()?;

    Ok(Solution {
        maneuvers,
        total_dv: resp.total_dv,
        iterations: resp.iterations,
        residual: resp.residual,
        lambda: resp.lambda.to_vec(),
    })
}

/// Parse a JSON `SolveRequest`, run it, and return the JSON `SolveResponse`.
///
/// Raises `ValueError` for malformed JSON or invalid input, `RuntimeError` for
/// solver failures.
#[pyfunction]
fn solve_json(input: &str) -> PyResult<String> {
    koenig_planner_api::run_json(input).map_err(api_err_to_py)
}

/// The `koenig_planner` Python module.
#[pymodule]
fn koenig_planner(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<Orbit>()?;
    m.add_class::<Maneuver>()?;
    m.add_class::<Solution>()?;
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add_function(wrap_pyfunction!(solve_json, m)?)?;
    Ok(())
}
