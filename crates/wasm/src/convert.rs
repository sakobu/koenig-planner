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
            w_meters,
            cost,
            params,
            initial_times,
        } = r;
        api::SolveRequest {
            chief: chief.into(),
            t_i: *t_i,
            t_f: *t_f,
            dt: *dt,
            w_meters: *w_meters,
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
        let api::SolveResponse {
            maneuvers,
            total_dv,
            iterations,
            residual,
            lambda,
            primer_times,
            primer_magnitude,
            primer_rtn,
        } = resp;
        dto::SolveResponse {
            maneuvers: maneuvers.iter().map(Into::into).collect(),
            total_dv,
            iterations,
            residual,
            lambda,
            primer_times,
            primer_magnitude,
            primer_rtn,
            geometry, // presentation-only field merged in
        }
    }
}

impl From<api::ApiErrorKind> for dto::ApiErrorKind {
    fn from(k: api::ApiErrorKind) -> Self {
        match k {
            api::ApiErrorKind::BadRequest => dto::ApiErrorKind::BadRequest,
            api::ApiErrorKind::Solver => dto::ApiErrorKind::Solver,
            api::ApiErrorKind::Internal => dto::ApiErrorKind::Internal,
        }
    }
}

impl From<api::ApiError> for dto::ApiError {
    fn from(e: api::ApiError) -> Self {
        dto::ApiError {
            kind: e.kind.into(),
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
            w_meters: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
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
            w_meters: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
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
            kind: api::ApiErrorKind::BadRequest,
            message: "boom".to_string(),
        };
        let m: dto::ApiError = e.into();
        assert_eq!(m.kind, dto::ApiErrorKind::BadRequest);
        assert_eq!(m.message, "boom");
    }

    #[wasm_bindgen_test]
    fn response_converts_field_for_field() {
        let resp = api::SolveResponse {
            maneuvers: vec![api::ManeuverDto {
                t: 7.0,
                dv: [1.0, 2.0, 3.0],
            }],
            total_dv: 4.0,
            iterations: 2,
            residual: 1e-9,
            lambda: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            primer_times: vec![0.0, 7.0],
            primer_magnitude: vec![0.4, 1.0],
            primer_rtn: vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
        };
        let geom = dto::ChiefGeometry {
            a: 25_000e3,
            e: 0.7,
            maneuver_nu: vec![0.5],
            perigee_window: Some([0.1, 0.2]),
            orbit_eci: vec![[1.0, 2.0, 3.0]],
            chief_track_eci: vec![[4.0, 5.0, 6.0]],
            maneuver_eci: vec![dto::ManeuverEciDto {
                position_eci: [7.0, 8.0, 9.0],
                dv_eci: [0.1, 0.2, 0.3],
            }],
            maneuver_rtn: vec![dto::ManeuverRtnDto {
                position_rtn: [0.7, 0.8, 0.9],
                dv_rtn: [0.01, 0.02, 0.03],
            }],
            primer_eci: vec![[0.4, 0.5, 0.6]],
            primer_rtn: vec![[0.7, 0.8, 0.9]],
            perigee_arc_eci: Some(vec![[1.1, 1.2, 1.3]]),
            deputy_track_rtn: vec![[1.0, 2.0, 3.0]],
            target_roe: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
        };
        let got: dto::SolveResponse = (resp, geom).into();
        assert_eq!(got.maneuvers.len(), 1);
        assert_eq!(got.maneuvers[0].t, 7.0);
        assert_eq!(got.maneuvers[0].dv, [1.0, 2.0, 3.0]);
        assert_eq!(got.total_dv, 4.0);
        assert_eq!(got.iterations, 2);
        assert_eq!(got.residual, 1e-9);
        assert_eq!(got.lambda, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(got.primer_times, vec![0.0, 7.0]);
        assert_eq!(got.primer_magnitude, vec![0.4, 1.0]);
        assert_eq!(got.primer_rtn, vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]);
        assert_eq!(got.geometry.a, 25_000e3);
        assert_eq!(got.geometry.e, 0.7);
        assert_eq!(got.geometry.maneuver_nu, vec![0.5]);
        assert_eq!(got.geometry.perigee_window, Some([0.1, 0.2]));
        assert_eq!(got.geometry.maneuver_rtn.len(), 1);
        assert_eq!(got.geometry.maneuver_rtn[0].position_rtn, [0.7, 0.8, 0.9]);
        assert_eq!(got.geometry.maneuver_rtn[0].dv_rtn, [0.01, 0.02, 0.03]);
        assert_eq!(got.geometry.primer_rtn, vec![[0.7, 0.8, 0.9]]);
    }
}
