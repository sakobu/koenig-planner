//! Conversions between the wasm mirror DTOs and the `crates/api` DTOs.
//!
//! Every conversion — request and response — destructures all fields (no `..`),
//! so adding or renaming a field in `crates/api` breaks compilation here until
//! the mirror is updated: drift is impossible to merge in either direction. The
//! response conversion takes `(api::SolveResponse, ChiefGeometry)` as a tuple
//! because the api response has no presentation `geometry` field.

use crate::dto;
use koenig_damico_planner_api as api;

impl From<&dto::OrbitDto> for api::OrbitDto {
    fn from(o: &dto::OrbitDto) -> Self {
        let dto::OrbitDto {
            a,
            e,
            i,
            raan,
            argp,
            mean_anom,
        } = *o;
        api::OrbitDto {
            a,
            e,
            i,
            raan,
            argp,
            mean_anom,
        }
    }
}

impl From<&dto::CostSpec> for api::CostSpec {
    fn from(c: &dto::CostSpec) -> Self {
        match c {
            dto::CostSpec::Norm2 => api::CostSpec::Norm2,
            dto::CostSpec::FaceMax => api::CostSpec::FaceMax,
            dto::CostSpec::Piecewise { period, t_perigee0 } => api::CostSpec::Piecewise {
                period: *period,
                t_perigee0: *t_perigee0,
            },
        }
    }
}

impl From<&dto::SolveParamsDto> for api::SolveParamsDto {
    fn from(p: &dto::SolveParamsDto) -> Self {
        let dto::SolveParamsDto {
            n_coarse,
            n_init,
            eps_cost,
            eps_remove,
        } = *p;
        api::SolveParamsDto {
            n_coarse,
            n_init,
            eps_cost,
            eps_remove,
        }
    }
}

impl From<&dto::SolveRequest> for api::SolveRequest {
    fn from(r: &dto::SolveRequest) -> Self {
        let dto::SolveRequest {
            chief,
            t_i,
            t_f,
            dt,
            w_metres,
            cost,
            params,
            initial_times,
        } = r;
        api::SolveRequest {
            chief: chief.into(),
            t_i: *t_i,
            t_f: *t_f,
            dt: *dt,
            w_metres: *w_metres,
            cost: cost.into(),
            params: params.as_ref().map(Into::into),
            initial_times: initial_times.clone(),
        }
    }
}

impl From<&api::ManeuverDto> for dto::ManeuverDto {
    fn from(m: &api::ManeuverDto) -> Self {
        let api::ManeuverDto { t, dv } = m; // api::ManeuverDto is not Copy
        dto::ManeuverDto { t: *t, dv: *dv }
    }
}

impl From<(api::SolveResponse, dto::ChiefGeometry)> for dto::SolveResponse {
    fn from((resp, geometry): (api::SolveResponse, dto::ChiefGeometry)) -> Self {
        let api::SolveResponse { maneuvers, total_dv, iterations, residual, lambda } = resp;
        dto::SolveResponse {
            maneuvers: maneuvers.iter().map(Into::into).collect(),
            total_dv,
            iterations,
            residual,
            lambda,
            geometry, // presentation-only field merged in
        }
    }
}

impl From<api::ApiError> for dto::ApiError {
    fn from(e: api::ApiError) -> Self {
        dto::ApiError {
            kind: e.kind.to_string(),
            message: e.message,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wasm_bindgen_test::wasm_bindgen_test;

    fn golden() -> dto::SolveRequest {
        dto::SolveRequest {
            chief: dto::OrbitDto {
                a: 25_000e3,
                e: 0.7,
                i: 40.0,
                raan: 358.0,
                argp: 0.0,
                mean_anom: 180.0,
            },
            t_i: 0.0,
            t_f: 117_990.0,
            dt: 30.0,
            w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
            cost: dto::CostSpec::Piecewise {
                period: None,
                t_perigee0: None,
            },
            params: None,
            initial_times: None,
        }
    }

    #[wasm_bindgen_test]
    fn request_converts_field_for_field() {
        let got: api::SolveRequest = (&golden()).into();
        let want = api::SolveRequest {
            chief: api::OrbitDto {
                a: 25_000e3,
                e: 0.7,
                i: 40.0,
                raan: 358.0,
                argp: 0.0,
                mean_anom: 180.0,
            },
            t_i: 0.0,
            t_f: 117_990.0,
            dt: 30.0,
            w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
            cost: api::CostSpec::Piecewise {
                period: None,
                t_perigee0: None,
            },
            params: None,
            initial_times: None,
        };
        assert_eq!(got, want);
    }

    #[wasm_bindgen_test]
    fn api_error_maps_kind_and_message() {
        let e = api::ApiError {
            kind: "bad_request",
            message: "boom".to_string(),
        };
        let m: dto::ApiError = e.into();
        assert_eq!(m.kind, "bad_request");
        assert_eq!(m.message, "boom");
    }

    #[wasm_bindgen_test]
    fn response_converts_field_for_field() {
        let resp = api::SolveResponse {
            maneuvers: vec![api::ManeuverDto { t: 7.0, dv: [1.0, 2.0, 3.0] }],
            total_dv: 4.0,
            iterations: 2,
            residual: 1e-9,
            lambda: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        };
        let geom = dto::ChiefGeometry {
            a: 25_000e3,
            e: 0.7,
            maneuver_nu: vec![0.5],
            perigee_window: Some([0.1, 0.2]),
        };
        let got: dto::SolveResponse = (resp, geom).into();
        assert_eq!(got.maneuvers.len(), 1);
        assert_eq!(got.maneuvers[0].t, 7.0);
        assert_eq!(got.maneuvers[0].dv, [1.0, 2.0, 3.0]);
        assert_eq!(got.total_dv, 4.0);
        assert_eq!(got.iterations, 2);
        assert_eq!(got.residual, 1e-9);
        assert_eq!(got.lambda, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(got.geometry.a, 25_000e3);
        assert_eq!(got.geometry.maneuver_nu, vec![0.5]);
        assert_eq!(got.geometry.perigee_window, Some([0.1, 0.2]));
    }
}
