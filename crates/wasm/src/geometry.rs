//! Presentation geometry for the orbit panel, computed by REUSING the core's
//! FD-verified Kepler solver and J2-secular propagation — no math is
//! re-implemented here.

use crate::dto;
use crate::frames;
use koenig_damico_planner_api as api;
use koenig_damico_planner_api::core::cost::Piecewise;
use koenig_damico_planner_api::core::dynamics::AbsoluteOrbit;
use koenig_damico_planner_api::core::PlannerError;
use std::f64::consts::TAU;

/// Sample counts for the presentation curves (chief orbit loop, perigee arc).
const N_ORBIT_SAMPLES: usize = 256;
const N_ARC_SAMPLES: usize = 64;

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

/// Presentation geometry for the 3D scene + orbit panel, computed by REUSING
/// the core's Kepler solver and J2-secular propagation. Reads maneuver times,
/// Δv, and the primer history from the api response. The deputy's relative
/// motion (RTN frame, over the playback grid) is reconstructed via exact ROE
/// inversion and propagated with the same core dynamics as the chief.
///
/// # Errors
/// Propagates the core `true_anomaly`/`Piecewise` errors (non-elliptic `e`).
pub fn chief_geometry(
    req: &dto::SolveRequest,
    resp: &api::SolveResponse,
) -> Result<dto::ChiefGeometry, PlannerError> {
    let chief = chief_orbit(&req.chief);
    // Maneuver and primer grid times are ABSOLUTE (`t_i + k·dt`, [KD20]); the core's
    // `propagate` advances the `t_i` epoch by a DURATION, so evaluate at `t - t_i`.
    let dur = |t: f64| t - req.t_i;

    // True anomaly at each maneuver time.
    let mut maneuver_nu = Vec::with_capacity(resp.maneuvers.len());
    for m in &resp.maneuvers {
        maneuver_nu.push(chief.propagate(dur(m.t)).true_anomaly()?);
    }

    // Closed-loop chief-orbit shape in ECI (sampled by evenly-spaced true
    // anomaly using the chief's epoch elements; the orbit precesses only
    // slowly over the window).
    let mut orbit_eci = Vec::with_capacity(N_ORBIT_SAMPLES + 1);
    for k in 0..=N_ORBIT_SAMPLES {
        let nu = TAU * (k as f64) / (N_ORBIT_SAMPLES as f64);
        orbit_eci.push(frames::orbit_point_eci(
            chief.a, chief.e, chief.i, chief.raan, chief.argp, nu,
        ));
    }

    // Chief position at each primer sample (playback track).
    let mut chief_track_eci = Vec::with_capacity(resp.primer_times.len());
    for &t in &resp.primer_times {
        chief_track_eci.push(frames::position_eci(&chief.propagate(dur(t)))?);
    }

    // Burn position + Δv direction in ECI.
    let mut maneuver_eci = Vec::with_capacity(resp.maneuvers.len());
    for m in &resp.maneuvers {
        let orb = chief.propagate(dur(m.t));
        maneuver_eci.push(dto::ManeuverEciDto {
            position_eci: frames::position_eci(&orb)?,
            dv_eci: frames::rtn_to_eci(&orb, m.dv)?,
        });
    }

    // Primer vector in ECI at each primer sample (RTN→ECI at the chief there).
    let mut primer_eci = Vec::with_capacity(resp.primer_times.len());
    for (&t, &p_rtn) in resp.primer_times.iter().zip(resp.primer_rtn.iter()) {
        primer_eci.push(frames::rtn_to_eci(&chief.propagate(dur(t)), p_rtn)?);
    }

    let perigee_window = match &req.cost {
        dto::CostSpec::Piecewise { period, t_perigee0 } => {
            let period = period.unwrap_or(TAU / chief.mean_motion());
            // Duration from the t_i epoch to perigee; the same default the solver
            // uses (api::run), so the drawn window matches the applied cost.
            let t_pc = t_perigee0.unwrap_or_else(|| chief.time_to_perigee());
            let pw = match t_perigee0 {
                Some(tp) => Piecewise::with_perigee_epoch(period, *tp),
                None => Piecewise::with_perigee_epoch(period, t_pc),
            }?;
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

    // ECI samples of the perigee-window arc (piecewise only).
    let perigee_arc_eci = perigee_window.map(|[lo, hi]| {
        (0..=N_ARC_SAMPLES)
            .map(|k| {
                let nu = lo + (hi - lo) * (k as f64) / (N_ARC_SAMPLES as f64);
                frames::orbit_point_eci(chief.a, chief.e, chief.i, chief.raan, chief.argp, nu)
            })
            .collect()
    });

    // Deputy relative orbit: reconstruct the deputy's absolute orbit from the
    // (dimensionless) target ROE via the exact ROE inverse. It is propagated
    // with the core below — at each playback sample (deputy_track_rtn) and at
    // each maneuver time — then differenced in ECI and expressed in the chief
    // RTN frame. Faithful by reuse — only the ROE inverse and frame rotations
    // are new.
    let roe = [
        req.w_meters[0] / chief.a,
        req.w_meters[1] / chief.a,
        req.w_meters[2] / chief.a,
        req.w_meters[3] / chief.a,
        req.w_meters[4] / chief.a,
        req.w_meters[5] / chief.a,
    ];
    let deputy = frames::deputy_from_roe(&chief, roe);

    // Deputy position in chief RTN at each playback (primer_times) sample, so the
    // RTN view's glyph tracks the same scrubber as the ECI spacecraft.
    let mut deputy_track_rtn = Vec::with_capacity(resp.primer_times.len());
    for &t in &resp.primer_times {
        let c_t = chief.propagate(dur(t));
        let d_t = deputy.propagate(dur(t));
        let r_c = frames::position_eci(&c_t)?;
        let r_d = frames::position_eci(&d_t)?;
        let rel_eci = [r_d[0] - r_c[0], r_d[1] - r_c[1], r_d[2] - r_c[2]];
        deputy_track_rtn.push(frames::eci_to_rtn(&c_t, rel_eci)?);
    }

    // Burn position + native-RTN Δv per maneuver, in the chief RTN frame (the
    // RTN analog of maneuver_eci). Anchor only: position_rtn is the deputy's
    // relative RTN location at the burn time, reconstructed exactly like
    // deputy_track_rtn (same target-ROE deputy), so each marker sits on the
    // drawn deputy curve. The real initial→target transition arc isn't
    // reconstructible from (t, Δv), so this position is schematic and can differ
    // visibly at meters scale; only the Δv DIRECTION is exact (magnitude lives
    // in the R/T/N Δv-component bars). dv_rtn is m.dv echoed with no rotation —
    // it is already in the chief RTN frame (like target_roe echoes w_meters).
    let mut maneuver_rtn = Vec::with_capacity(resp.maneuvers.len());
    for m in &resp.maneuvers {
        let c_t = chief.propagate(dur(m.t));
        let d_t = deputy.propagate(dur(m.t));
        let r_c = frames::position_eci(&c_t)?;
        let r_d = frames::position_eci(&d_t)?;
        let rel_eci = [r_d[0] - r_c[0], r_d[1] - r_c[1], r_d[2] - r_c[2]];
        maneuver_rtn.push(dto::ManeuverRtnDto {
            position_rtn: frames::eci_to_rtn(&c_t, rel_eci)?,
            dv_rtn: m.dv,
        });
    }

    // Primer vector in the chief RTN frame at each playback sample — presentation
    // copy of resp.primer_rtn (mirrors primer_eci) so the RTN scene draws the
    // swept primer arrow from geometry alone.
    let primer_rtn = resp.primer_rtn.clone();

    Ok(dto::ChiefGeometry {
        a: req.chief.a,
        e: req.chief.e,
        maneuver_nu,
        perigee_window,
        orbit_eci,
        chief_track_eci,
        maneuver_eci,
        maneuver_rtn,
        primer_eci,
        primer_rtn,
        perigee_arc_eci,
        deputy_track_rtn,
        target_roe: req.w_meters,
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
            w_meters: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
            cost,
            params: None,
            initial_times: None,
        }
    }

    // Minimal api response carrying the maneuver/primer data chief_geometry reads.
    fn resp_with(maneuver_times: &[f64]) -> api::SolveResponse {
        api::SolveResponse {
            maneuvers: maneuver_times
                .iter()
                .map(|&t| api::ManeuverDto {
                    t,
                    dv: [1.0, 0.0, 0.0],
                })
                .collect(),
            total_dv: 0.0,
            iterations: 0,
            residual: 0.0,
            lambda: [0.0; 6],
            primer_times: vec![0.0, 1000.0, 2000.0],
            primer_magnitude: vec![0.5, 1.0, 0.5],
            primer_rtn: vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        }
    }

    #[wasm_bindgen_test]
    fn perigee_at_epoch_gives_zero_true_anomaly() {
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp_with(&[0.0])).unwrap();
        assert_eq!(g.maneuver_nu.len(), 1);
        assert!(
            g.maneuver_nu[0].abs() < 1e-9,
            "ν at perigee ~0, got {}",
            g.maneuver_nu[0]
        );
        assert!(g.perigee_window.is_none(), "norm2 has no perigee window");
        assert!(g.perigee_arc_eci.is_none(), "norm2 has no perigee arc");
    }

    #[wasm_bindgen_test]
    fn geometry_uses_relative_epoch_for_nonzero_t_i() {
        // Chief at perigee at its t_i epoch (M0 = 0). A burn exactly at t = t_i is a
        // zero-duration propagation, so ν ≈ 0 — independently of the absolute t_i.
        // With the pre-fix absolute-time propagation this lands ~1.7 rad off perigee.
        let mut req = req_with(dto::CostSpec::Norm2, 0.0);
        req.t_i = 50_000.0;
        req.t_f = 167_990.0;
        let mut resp = resp_with(&[req.t_i]); // one maneuver at t = t_i
        resp.primer_times = vec![req.t_i, req.t_i + 1000.0];
        resp.primer_rtn = vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let g = chief_geometry(&req, &resp).unwrap();
        assert!(
            g.maneuver_nu[0].abs() < 1e-9,
            "burn at t=t_i must be ν≈0 (relative epoch), got {}",
            g.maneuver_nu[0]
        );
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
            &resp_with(&[]),
        )
        .unwrap();
        let [lo, hi] = g.perigee_window.expect("piecewise has a window");
        assert!(
            lo < 0.0 && hi > 0.0,
            "window should bracket perigee, got [{lo}, {hi}]"
        );
        assert!(lo > -PI && hi < PI);
        assert!(g.perigee_arc_eci.is_some(), "piecewise has a perigee arc");
    }

    #[wasm_bindgen_test]
    fn piecewise_window_brackets_perigee_general_m0() {
        let g = chief_geometry(
            &req_with(
                dto::CostSpec::Piecewise {
                    period: None,
                    t_perigee0: None,
                },
                90.0,
            ),
            &resp_with(&[]),
        )
        .unwrap();
        let [lo, hi] = g.perigee_window.expect("piecewise has a window");
        assert!(
            lo < 0.0 && hi > 0.0,
            "window should bracket perigee for M0=90, got [{lo}, {hi}]"
        );
        assert!(lo > -PI && hi < PI);
    }

    #[wasm_bindgen_test]
    fn orbit_loop_is_closed_and_radius_bounded() {
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp_with(&[0.0])).unwrap();
        assert_eq!(g.orbit_eci.len(), N_ORBIT_SAMPLES + 1);
        let first = g.orbit_eci[0];
        let last = g.orbit_eci[N_ORBIT_SAMPLES];
        for k in 0..3 {
            assert!((first[k] - last[k]).abs() < 1.0, "loop not closed at k={k}");
        }
        let (a, e) = (25_000e3_f64, 0.7_f64);
        for p in &g.orbit_eci {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(
                r >= a * (1.0 - e) - 1.0 && r <= a * (1.0 + e) + 1.0,
                "r={r}"
            );
        }
    }

    #[wasm_bindgen_test]
    fn maneuver_and_primer_eci_preserve_magnitude_and_length() {
        let g = chief_geometry(
            &req_with(dto::CostSpec::Norm2, 0.0),
            &resp_with(&[0.0, 5000.0]),
        )
        .unwrap();
        // two maneuvers, each dv = [1,0,0] (RTN) → |dv_eci| == 1.
        assert_eq!(g.maneuver_eci.len(), 2);
        for m in &g.maneuver_eci {
            let mag = frames_norm(m.dv_eci);
            assert!((mag - 1.0).abs() < 1e-9, "|dv_eci| should be 1, got {mag}");
        }
        // primer_eci parallels primer_times (len 3); each primer_rtn = [0,1,0] → |p|=1.
        assert_eq!(g.primer_eci.len(), 3);
        assert_eq!(g.chief_track_eci.len(), 3);
        for p in &g.primer_eci {
            assert!((frames_norm(*p) - 1.0).abs() < 1e-9);
        }
    }

    #[wasm_bindgen_test]
    fn target_roe_echoes_w_meters() {
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp_with(&[])).unwrap();
        assert_eq!(g.target_roe, [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]);
    }

    #[wasm_bindgen_test]
    fn deputy_track_rtn_length_matches_primer_times() {
        // resp_with gives primer_times of length 3.
        let resp = resp_with(&[0.0]);
        assert_eq!(resp.primer_times.len(), 3);
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp).unwrap();
        assert_eq!(
            g.deputy_track_rtn.len(),
            resp.primer_times.len(),
            "deputy_track_rtn must be parallel to primer_times"
        );
    }

    #[wasm_bindgen_test]
    fn zero_roe_gives_zero_deputy_track_rtn() {
        let mut req = req_with(dto::CostSpec::Norm2, 0.0);
        req.w_meters = [0.0; 6];
        let g = chief_geometry(&req, &resp_with(&[])).unwrap();
        for p in &g.deputy_track_rtn {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(r < 1e-3, "zero ROE ⇒ coincident orbits in track, got r={r}");
        }
    }

    #[wasm_bindgen_test]
    fn maneuver_rtn_parallels_maneuvers_and_echoes_native_dv() {
        let g = chief_geometry(
            &req_with(dto::CostSpec::Norm2, 0.0),
            &resp_with(&[0.0, 5000.0]),
        )
        .unwrap();
        // One RTN marker per maneuver, matching the ECI markers.
        assert_eq!(g.maneuver_rtn.len(), 2);
        assert_eq!(g.maneuver_rtn.len(), g.maneuver_eci.len());
        // dv_rtn is m.dv echoed with NO rotation (resp_with uses dv = [1,0,0]).
        for m in &g.maneuver_rtn {
            assert_eq!(m.dv_rtn, [1.0, 0.0, 0.0]);
        }
    }

    #[wasm_bindgen_test]
    fn maneuver_rtn_position_is_meters_scale() {
        // Meters-scale ROE ⇒ the relative burn position is meters-scale, finite,
        // and far below the ~2.5e7 m chief radius (same bound as the rel orbit).
        let g = chief_geometry(
            &req_with(dto::CostSpec::Norm2, 0.0),
            &resp_with(&[0.0, 5000.0]),
        )
        .unwrap();
        for m in &g.maneuver_rtn {
            let p = m.position_rtn;
            assert!(p[0].is_finite() && p[1].is_finite() && p[2].is_finite());
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(r > 0.0 && r < 2.5e5, "relative burn scale off: {r}");
        }
    }

    #[wasm_bindgen_test]
    fn zero_roe_gives_zero_maneuver_rtn_position() {
        let mut req = req_with(dto::CostSpec::Norm2, 0.0);
        req.w_meters = [0.0; 6];
        let g = chief_geometry(&req, &resp_with(&[0.0, 5000.0])).unwrap();
        for m in &g.maneuver_rtn {
            let p = m.position_rtn;
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(r < 1e-3, "zero ROE ⇒ coincident orbits, got r={r}");
        }
    }

    #[wasm_bindgen_test]
    fn primer_rtn_echoes_response_at_each_sample() {
        let resp = resp_with(&[0.0]);
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp).unwrap();
        // Presentation copy: parallel to primer_times and byte-equal to the
        // response primer_rtn (RTN analog of primer_eci; no rotation).
        assert_eq!(g.primer_rtn.len(), resp.primer_times.len());
        assert_eq!(g.primer_rtn, resp.primer_rtn);
    }

    // Local norm helper (frames::norm is private to its module).
    fn frames_norm(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }
}
