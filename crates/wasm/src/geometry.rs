//! Presentation geometry for the orbit panel, computed by REUSING the core's
//! FD-verified Kepler solver and J2-secular propagation — no math is
//! re-implemented here.

use crate::dto;
use crate::frames;
use crate::roe_track;
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

/// Deputy position relative to the chief in the chief RTN frame, both orbits
/// already evaluated at the same instant (callers propagate each from its own
/// epoch).
///
/// # Errors
/// Fails when either orbit is non-elliptic (`e ≥ 1`), which has no Kepler
/// solution. Callers treat this as "relative track not drawable" (best-effort) —
/// it is a presentation artifact, not a solver failure, so it must not sink an
/// otherwise-valid solve.
fn rel_rtn(c_t: &AbsoluteOrbit, d_t: &AbsoluteOrbit) -> Result<[f64; 3], PlannerError> {
    let r_c = frames::position_eci(c_t)?;
    let r_d = frames::position_eci(d_t)?;
    let rel_eci = [r_d[0] - r_c[0], r_d[1] - r_c[1], r_d[2] - r_c[2]];
    frames::eci_to_rtn(c_t, rel_eci)
}

/// Index of the `times` sample nearest `t` (monotonic grid; burn times are grid
/// members, so this is exact up to float drift).
fn nearest_sample(times: &[f64], t: f64) -> usize {
    let mut best = 0;
    let mut best_d = f64::INFINITY;
    for (k, &tk) in times.iter().enumerate() {
        let d = (tk - t).abs();
        if d < best_d {
            best_d = d;
            best = k;
        }
    }
    best
}

/// Presentation geometry for the 3D scene + orbit panel, computed by REUSING
/// the core's Kepler solver and J2-secular propagation. Reads maneuver times,
/// Δv, and the primer history from the api response. The deputy's relative
/// motion (RTN frame, over the playback grid) is reconstructed via exact ROE
/// inversion and propagated with the same core dynamics as the chief.
///
/// # Errors
/// Propagates the core `true_anomaly` / `Piecewise` errors for the **chief**
/// (e.g. a non-elliptic chief) — unreachable once the solver has accepted the
/// request. The **deputy**-derived fields (`target_track_rtn`,
/// `transfer_track_rtn`, `maneuver_rtn`) are best-effort: they degrade to
/// empty for a non-elliptic reconstructed deputy rather than failing the
/// whole geometry (see [`rel_rtn`]).
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

    // Chief true anomaly at each playback sample (scrub readout).
    let mut chief_nu_track = Vec::with_capacity(resp.primer_times.len());
    for &t in &resp.primer_times {
        chief_nu_track.push(chief.propagate(dur(t)).true_anomaly()?);
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

    // Target deputy: the orbit whose ROE relative to the chief AT t_f is the
    // request target — the epoch where the solver enforces δα(t_f) = w/a —
    // reconstructed via the exact ROE inverse. Samples before t_f propagate it
    // backward (negative durations are exact: linear secular rates + Kepler).
    let roe = [
        req.w_meters[0] / chief.a,
        req.w_meters[1] / chief.a,
        req.w_meters[2] / chief.a,
        req.w_meters[3] / chief.a,
        req.w_meters[4] / chief.a,
        req.w_meters[5] / chief.a,
    ];
    let deputy_tgt = frames::deputy_from_roe(&chief.propagate(req.t_f - req.t_i), roe);

    // Target relative orbit at each playback sample — the ghost curve the
    // transfer lands on. Best-effort: empty when the target ROE implies a
    // non-elliptic deputy, so an extreme target degrades the relative track,
    // not the whole solve.
    let target_track_rtn: Vec<[f64; 3]> = resp
        .primer_times
        .iter()
        .map(|&t| rel_rtn(&chief.propagate(dur(t)), &deputy_tgt.propagate(t - req.t_f)))
        .collect::<Result<_, _>>()
        .unwrap_or_default();

    // Primer vector in the chief RTN frame at each playback sample — presentation
    // copy of resp.primer_rtn (mirrors primer_eci) so the RTN scene draws the
    // swept primer arrow from geometry alone.
    let primer_rtn = resp.primer_rtn.clone();

    // Controlled mean-ROE trajectory + per-burn jumps ([KD20] eq. 11): the
    // pseudostate accumulated from δα = 0 at t_i by the plan's burns, on the
    // playback grid. Chief-only reuse of the core's Φ/B — independent of the
    // reconstructed deputy, so it stays available when the deputy-derived
    // fields above degrade. Best-effort for idiom uniformity only (a chief
    // that already solved cannot fail B here).
    let (roe_track, roe_jumps) =
        roe_track::controlled_roe_track(&chief, req.t_i, &resp.maneuvers, &resp.primer_times)
            .unwrap_or_default();

    // True transfer trajectory: the controlled pseudostate δα(t) (roe_track, in
    // meters) mapped through the exact ROE inverse at each sample's chief and
    // differenced in ECI — the solver's own model drawn in position space. No
    // deputy propagation: δα(t) is already the state AT t. Best-effort like the
    // target track (an extreme mid-transfer state can be non-elliptic); also
    // empty whenever roe_track itself degraded.
    let transfer_track_rtn: Vec<[f64; 3]> = resp
        .primer_times
        .iter()
        .zip(&roe_track)
        .map(|(&t, roe_m)| {
            let c_t = chief.propagate(dur(t));
            let d_t = frames::deputy_from_roe(&c_t, roe_m.map(|v| v / chief.a));
            rel_rtn(&c_t, &d_t)
        })
        .collect::<Result<_, _>>()
        .unwrap_or_default();

    // Burn markers ride the transfer: each position is the transfer sample at
    // the burn time (post-burn, roe_track's inclusive convention), so markers
    // sit bitwise on the drawn polyline's kinks. Burn times are grid members
    // (the solver extracts them from the grid); nearest-sample lookup absorbs
    // float drift. dv_rtn is m.dv echoed — already the chief RTN frame. Empty
    // whenever the transfer itself is unavailable.
    let maneuver_rtn: Vec<dto::ManeuverRtnDto> = if !resp.primer_times.is_empty()
        && transfer_track_rtn.len() == resp.primer_times.len()
    {
        resp.maneuvers
            .iter()
            .map(|m| dto::ManeuverRtnDto {
                position_rtn: transfer_track_rtn[nearest_sample(&resp.primer_times, m.t)],
                dv_rtn: m.dv,
            })
            .collect()
    } else {
        Vec::new()
    };

    Ok(dto::ChiefGeometry {
        a: req.chief.a,
        e: req.chief.e,
        maneuver_nu,
        chief_nu_track,
        perigee_window,
        orbit_eci,
        chief_track_eci,
        maneuver_eci,
        maneuver_rtn,
        primer_eci,
        primer_rtn,
        perigee_arc_eci,
        target_track_rtn,
        transfer_track_rtn,
        roe_track,
        roe_jumps,
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
    fn target_track_rtn_length_matches_primer_times() {
        // resp_with gives primer_times of length 3.
        let resp = resp_with(&[0.0]);
        assert_eq!(resp.primer_times.len(), 3);
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp).unwrap();
        assert_eq!(
            g.target_track_rtn.len(),
            resp.primer_times.len(),
            "target_track_rtn must be parallel to primer_times"
        );
    }

    #[wasm_bindgen_test]
    fn zero_roe_gives_zero_target_track_rtn() {
        let mut req = req_with(dto::CostSpec::Norm2, 0.0);
        req.w_meters = [0.0; 6];
        let g = chief_geometry(&req, &resp_with(&[])).unwrap();
        for p in &g.target_track_rtn {
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

    #[wasm_bindgen_test]
    fn negative_duration_propagate_inverts_forward() {
        // Sampling the t_f-anchored target before t_f propagates backward:
        // a, e, i are constant and Ω/ω/M are linear at fixed secular rates, so
        // propagate(dt) then propagate(−dt) must return to the epoch elements.
        let c = chief_orbit(&req_with(dto::CostSpec::Norm2, 0.0).chief);
        let back = c.propagate(5000.0).propagate(-5000.0);
        assert!((back.a - c.a).abs() < 1e-6);
        assert!((back.raan - c.raan).abs() < 1e-12);
        assert!((back.argp - c.argp).abs() < 1e-12);
        assert!((back.mean_anom - c.mean_anom).abs() < 1e-12);
    }

    #[wasm_bindgen_test]
    fn target_track_is_anchored_at_t_f() {
        // δa ≠ 0 makes the anchor epoch observable (≈1.4 km of along-track drift
        // over this window): the last grid sample must sit exactly at the
        // t_f-anchored reconstruction, not the t_i-anchored one.
        let req = req_with(dto::CostSpec::Norm2, 0.0);
        let mut resp = resp_with(&[]);
        resp.primer_times = vec![0.0, 58_995.0, 117_990.0]; // ends exactly at t_f
        resp.primer_rtn = vec![[0.0, 1.0, 0.0]; 3];
        let g = chief_geometry(&req, &resp).unwrap();
        let chief_tf = chief_orbit(&req.chief).propagate(req.t_f - req.t_i);
        let roe = req.w_meters.map(|w| w / req.chief.a);
        let deputy_tf = frames::deputy_from_roe(&chief_tf, roe);
        let want = rel_rtn(&chief_tf, &deputy_tf).unwrap();
        let got = g.target_track_rtn[2];
        for k in 0..3 {
            assert!((got[k] - want[k]).abs() < 1e-6, "k={k}: got {} want {}", got[k], want[k]);
        }
    }

    #[wasm_bindgen_test]
    fn transfer_and_nu_tracks_parallel_primer_times() {
        let resp = resp_with(&[1000.0]);
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp).unwrap();
        assert_eq!(g.transfer_track_rtn.len(), resp.primer_times.len());
        assert_eq!(g.chief_nu_track.len(), resp.primer_times.len());
    }

    #[wasm_bindgen_test]
    fn transfer_is_zero_before_the_first_burn() {
        // δα = 0 exactly before the first burn, but the inversion round-trips
        // atan2/sqrt, so the reconstructed deputy matches the chief only to
        // fp round-off — sub-µm at a = 2.5e7 m, not bitwise zero.
        let resp = resp_with(&[2000.0]); // burn at the LAST grid sample
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp).unwrap();
        for k in 0..2 {
            let r = frames_norm(g.transfer_track_rtn[k]);
            assert!(r < 1e-6, "pre-burn sample {k} should be ~origin, got {r} m");
        }
    }

    #[wasm_bindgen_test]
    fn transfer_matches_the_independent_chain() {
        // Pins units (meters → ÷a), the per-sample chief epoch, and the
        // post-burn-inclusive convention against a from-parts recomputation.
        let req = req_with(dto::CostSpec::Norm2, 0.0);
        let resp = resp_with(&[1000.0]);
        let g = chief_geometry(&req, &resp).unwrap();
        let chief = chief_orbit(&req.chief);
        let (track_m, _) = roe_track::controlled_roe_track(
            &chief,
            req.t_i,
            &resp.maneuvers,
            &resp.primer_times,
        )
        .unwrap();
        let k = 2;
        let c_t = chief.propagate(resp.primer_times[k] - req.t_i);
        let d_t = frames::deputy_from_roe(&c_t, track_m[k].map(|v| v / chief.a));
        let want = rel_rtn(&c_t, &d_t).unwrap();
        for (j, w) in want.iter().enumerate() {
            assert!((g.transfer_track_rtn[k][j] - w).abs() < 1e-9, "j={j}");
        }
    }

    #[wasm_bindgen_test]
    fn maneuver_markers_sit_on_the_transfer() {
        // Marker positions are the transfer samples at the burn times — bitwise,
        // so they can never float off the drawn polyline.
        let resp = resp_with(&[1000.0, 2000.0]);
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp).unwrap();
        assert_eq!(g.maneuver_rtn.len(), 2);
        assert_eq!(g.maneuver_rtn[0].position_rtn, g.transfer_track_rtn[1]);
        assert_eq!(g.maneuver_rtn[1].position_rtn, g.transfer_track_rtn[2]);
    }

    #[wasm_bindgen_test]
    fn nu_track_zero_at_perigee_epoch() {
        // M₀ = 0 chief and the first sample at t = t_i ⇒ ν ≈ 0.
        let g = chief_geometry(&req_with(dto::CostSpec::Norm2, 0.0), &resp_with(&[])).unwrap();
        assert!(g.chief_nu_track[0].abs() < 1e-9, "got {}", g.chief_nu_track[0]);
    }

    #[wasm_bindgen_test]
    fn burn_position_discontinuity_is_second_order_in_dv() {
        // An impulse changes velocity, not position: the GVE kick B·Δv preserves
        // the reconstructed position to first order, so the pre→post gap must
        // shrink ~4× when Δv halves (O(|Δv|²) convergence) — the FD-style check
        // that the kick ↔ position mapping is physically consistent.
        let chief = chief_orbit(&req_with(dto::CostSpec::Norm2, 0.0).chief);
        let gap = |scale: f64| -> f64 {
            let ms = vec![api::ManeuverDto {
                t: 1000.0,
                dv: [0.5 * scale, 0.3 * scale, 0.2 * scale],
            }];
            let times = vec![0.0, 1000.0, 2000.0];
            let (track, jumps) =
                roe_track::controlled_roe_track(&chief, 0.0, &ms, &times).unwrap();
            let c_t = chief.propagate(1000.0);
            let post =
                rel_rtn(&c_t, &frames::deputy_from_roe(&c_t, track[1].map(|v| v / chief.a)))
                    .unwrap();
            let mut pre_roe = [0.0; 6];
            for i in 0..6 {
                pre_roe[i] = (track[1][i] - jumps[0][i]) / chief.a;
            }
            let pre = rel_rtn(&c_t, &frames::deputy_from_roe(&c_t, pre_roe)).unwrap();
            frames_norm([post[0] - pre[0], post[1] - pre[1], post[2] - pre[2]])
        };
        let g1 = gap(1.0);
        let g2 = gap(0.5);
        assert!(g1 > 1e-4, "gap should be measurable, got {g1} m");
        let ratio = g1 / g2;
        assert!((3.0..5.0).contains(&ratio), "expected ~4× shrink, got {ratio}");
    }
}
