//! Conversions between the wasm mirror DTOs and the `crates/api` DTOs.
//!
//! Every struct conversion destructures all fields (no `..`), so adding or
//! renaming a field in `crates/api` breaks compilation here until the mirror is
//! updated — drift is impossible to merge.

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
}
