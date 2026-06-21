//! Presentation geometry for the orbit panel, computed by REUSING the core's
//! FD-verified Kepler solver and J2-secular propagation — no math is
//! re-implemented here.

use crate::dto;
use koenig_damico_planner_api::core::cost::Piecewise;
use koenig_damico_planner_api::core::dynamics::AbsoluteOrbit;
use koenig_damico_planner_api::core::PlannerError;
use std::f64::consts::TAU;

/// Build the chief orbit (degrees → radians) the same way `api::run` does.
fn chief_orbit(c: &dto::OrbitDto) -> AbsoluteOrbit {
    AbsoluteOrbit::new(
        c.a,
        c.e,
        c.i.to_radians(),
        c.raan.to_radians(),
        c.argp.to_radians(),
        c.mean_anom.to_radians(),
    )
}

/// True anomaly at each maneuver time, plus the FaceMax perigee band for the
/// piecewise cost. Every angle comes from the core's `true_anomaly`
/// (`mean_to_true`) — faithful by reuse.
pub fn chief_geometry(
    req: &dto::SolveRequest,
    maneuver_times: &[f64],
) -> Result<dto::ChiefGeometry, PlannerError> {
    let chief = chief_orbit(&req.chief);

    let mut maneuver_nu = Vec::with_capacity(maneuver_times.len());
    for &t in maneuver_times {
        // propagate advances M at the J2 secular rate (matches the planner);
        // true_anomaly solves Kepler via the core's verified solver.
        maneuver_nu.push(chief.propagate(t).true_anomaly()?);
    }

    let perigee_window = match &req.cost {
        dto::CostSpec::Piecewise { period, t_perigee0 } => {
            let period = period.unwrap_or(TAU / chief.mean_motion());
            let t_pc = t_perigee0.unwrap_or(period / 2.0);
            let pw = match t_perigee0 {
                Some(tp) => Piecewise::with_perigee_epoch(period, *tp),
                None => Piecewise::new(period),
            }?;
            // Probe the ACTUAL eq.-49 selector outward from the perigee center to
            // find its time half-width — no hard-coded constant, no private fields.
            let step = (period / 720.0).max(1.0);
            let mut half = 0.0;
            while half < period / 2.0 && pw.in_perigee_window(t_pc + half + step) {
                half += step;
            }
            let nu_lo = chief.propagate(t_pc - half).true_anomaly()?;
            let nu_hi = chief.propagate(t_pc + half).true_anomaly()?;
            Some([nu_lo, nu_hi])
        }
        _ => None,
    };

    Ok(dto::ChiefGeometry {
        a: req.chief.a,
        e: req.chief.e,
        maneuver_nu,
        perigee_window,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use wasm_bindgen_test::wasm_bindgen_test;

    fn req_with(cost: dto::CostSpec, mean_anom: f64) -> dto::SolveRequest {
        dto::SolveRequest {
            chief: dto::OrbitDto {
                a: 25_000e3,
                e: 0.7,
                i: 40.0,
                raan: 358.0,
                argp: 0.0,
                mean_anom,
            },
            t_i: 0.0,
            t_f: 117_990.0,
            dt: 30.0,
            w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
            cost,
            params: None,
            initial_times: None,
        }
    }

    #[wasm_bindgen_test]
    fn perigee_at_epoch_gives_zero_true_anomaly() {
        // mean_anom = 0 is perigee → ν(t=0) ≈ 0.
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &[0.0]).unwrap();
        assert_eq!(g.maneuver_nu.len(), 1);
        assert!(
            g.maneuver_nu[0].abs() < 1e-9,
            "ν at perigee should be ~0, got {}",
            g.maneuver_nu[0]
        );
        assert!(g.perigee_window.is_none(), "norm2 has no perigee window");
    }

    #[wasm_bindgen_test]
    fn piecewise_window_brackets_perigee() {
        let g = chief_geometry(
            &req_with(
                dto::CostSpec::Piecewise {
                    period: None,
                    t_perigee0: None,
                },
                180.0,
            ),
            &[],
        )
        .unwrap();
        let [lo, hi] = g.perigee_window.expect("piecewise has a window");
        // The band straddles perigee (ν = 0): lo < 0 < hi, within (-π, π).
        assert!(
            lo < 0.0 && hi > 0.0,
            "window should bracket perigee, got [{lo}, {hi}]"
        );
        assert!(lo > -PI && hi < PI);
    }
}
