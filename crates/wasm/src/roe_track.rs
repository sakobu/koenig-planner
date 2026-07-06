//! Controlled mean-ROE trajectory of the impulsive plan: the pseudostate
//! `δα(t) = Σ_{t_j ≤ t} Φ(t_j→t)·B(t_j)·Δv_j` accumulated from `δα = 0` at
//! `t_i` (`[KD20]` eq. 11 impulsive state evolution; eq. 51 ROE ordering), and
//! the instantaneous per-burn jumps `B(t_j)·Δv_j`. Faithful by reuse: only the
//! core's FD-verified `state_transition` and `control_input_matrix` are
//! evaluated — per (sample, burn) pair, exactly like the solver's `Γ(t)` — and
//! no dynamics are re-implemented here.

use koenig_damico_planner_api as api;
use koenig_damico_planner_api::core::dynamics::b_matrix::control_input_matrix;
use koenig_damico_planner_api::core::dynamics::stm::state_transition;
use koenig_damico_planner_api::core::dynamics::AbsoluteOrbit;
use koenig_damico_planner_api::core::PlannerError;

/// `(roe_track, roe_jumps)` for the plan, both in **meters** — the
/// dimensionless ROE scaled by `chief_ti.a`, like the request's `w_meters`.
///
/// `roe_track[k]` is `a·δα(times[k])` with burns applied inclusively
/// (`t_j ≤ t`), so a burn's own grid sample carries the post-burn value and
/// `roe_track[k_j] − roe_jumps[j]` is the exact pre-burn state. `roe_jumps[j]`
/// is `a·B(t_j)·Δv_j`. Times are absolute grid times on the `t_i` axis;
/// `chief_ti` is the chief at the `t_i` epoch (propagation uses durations
/// `t − t_i`).
///
/// # Errors
/// Propagates `control_input_matrix`'s non-elliptic error — unreachable for a
/// chief that already produced a solve (the solver evaluated `B` at every burn
/// time); callers keep the best-effort idiom regardless.
#[allow(clippy::type_complexity)]
pub fn controlled_roe_track(
    chief_ti: &AbsoluteOrbit,
    t_i: f64,
    maneuvers: &[api::ManeuverDto],
    times: &[f64],
) -> Result<(Vec<[f64; 6]>, Vec<[f64; 6]>), PlannerError> {
    let a = chief_ti.a;
    // Per burn: the chief at the burn epoch and the dimensionless kick
    // B(t_j)·Δv_j. The same per-pair evaluation the solver's Γ(t) uses — no
    // semigroup shortcut (the secular-coefficient STM composes only
    // approximately; direct evaluation is the faithful form).
    let mut kicks: Vec<(f64, AbsoluteOrbit, [f64; 6])> = Vec::with_capacity(maneuvers.len());
    for m in maneuvers {
        let orb_j = chief_ti.propagate(m.t - t_i);
        let b = control_input_matrix(&orb_j)?;
        let mut kick = [0.0; 6];
        for (i, out) in kick.iter_mut().enumerate() {
            *out = b[(i, 0)] * m.dv[0] + b[(i, 1)] * m.dv[1] + b[(i, 2)] * m.dv[2];
        }
        kicks.push((m.t, orb_j, kick));
    }
    let roe_jumps = kicks
        .iter()
        .map(|(_, _, kick)| kick.map(|v| v * a))
        .collect();
    let roe_track = times
        .iter()
        .map(|&t| {
            let mut acc = [0.0; 6];
            for (t_j, orb_j, kick) in &kicks {
                if *t_j <= t {
                    let phi = state_transition(orb_j, t - t_j);
                    for (i, out) in acc.iter_mut().enumerate() {
                        for (k, kv) in kick.iter().enumerate() {
                            *out += phi[(i, k)] * kv;
                        }
                    }
                }
            }
            acc.map(|v| v * a)
        })
        .collect();
    Ok((roe_track, roe_jumps))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wasm_bindgen_test::wasm_bindgen_test;

    const A: f64 = 25_000e3;

    fn chief() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            A,
            0.7,
            40f64.to_radians(),
            358f64.to_radians(),
            0.0,
            180f64.to_radians(),
        )
    }

    fn grid(t_i: f64, n: usize, dt: f64) -> Vec<f64> {
        (0..n).map(|k| t_i + k as f64 * dt).collect()
    }

    fn burn(t: f64, dv: [f64; 3]) -> api::ManeuverDto {
        api::ManeuverDto { t, dv }
    }

    /// `|got − want| ≤ 1e-12·max(|want|, 1)` per component — scale-aware.
    fn assert_close(got: [f64; 6], want: [f64; 6], ctx: &str) {
        for i in 0..6 {
            let tol = 1e-12 * want[i].abs().max(1.0);
            assert!(
                (got[i] - want[i]).abs() <= tol,
                "{ctx}: component {i}: got {} want {}",
                got[i],
                want[i]
            );
        }
    }

    #[wasm_bindgen_test]
    fn outputs_parallel_their_inputs() {
        let times = grid(0.0, 101, 30.0);
        let ms = vec![
            burn(300.0, [0.01, 0.02, -0.005]),
            burn(1500.0, [-0.003, 0.004, 0.0]),
        ];
        let (track, jumps) = controlled_roe_track(&chief(), 0.0, &ms, &times).unwrap();
        assert_eq!(track.len(), times.len());
        assert_eq!(jumps.len(), ms.len());
    }

    #[wasm_bindgen_test]
    fn track_is_exactly_zero_before_the_first_burn() {
        let times = grid(0.0, 101, 30.0);
        let ms = vec![burn(1500.0, [0.01, 0.02, -0.005])];
        let (track, _) = controlled_roe_track(&chief(), 0.0, &ms, &times).unwrap();
        let mut checked = 0;
        for (k, &t) in times.iter().enumerate() {
            if t < 1500.0 {
                assert_eq!(track[k], [0.0; 6], "sample {k} precedes the burn");
                checked += 1;
            }
        }
        assert_eq!(checked, 50, "the fixture must actually cover the prefix");
    }

    #[wasm_bindgen_test]
    fn first_burn_sample_equals_its_jump() {
        // Phi(0) = I and no prior burns, so the post-burn sample IS the jump.
        let times = grid(0.0, 101, 30.0);
        let ms = vec![burn(1500.0, [0.01, 0.02, -0.005])];
        let (track, jumps) = controlled_roe_track(&chief(), 0.0, &ms, &times).unwrap();
        let k = times.iter().position(|&t| t == 1500.0).unwrap();
        assert_close(track[k], jumps[0], "post-burn sample vs jump");
    }

    #[wasm_bindgen_test]
    fn two_burn_track_is_the_sum_of_single_burn_tracks() {
        // The accumulation is linear in the burns; catches indexing bugs and
        // exercises coast propagation between burns.
        let times = grid(0.0, 101, 30.0);
        let (both, _) = controlled_roe_track(
            &chief(),
            0.0,
            &[
                burn(900.0, [0.01, 0.02, -0.005]),
                burn(2400.0, [-0.003, 0.004, 0.001]),
            ],
            &times,
        )
        .unwrap();
        let (only1, _) =
            controlled_roe_track(&chief(), 0.0, &[burn(900.0, [0.01, 0.02, -0.005])], &times)
                .unwrap();
        let (only2, _) = controlled_roe_track(
            &chief(),
            0.0,
            &[burn(2400.0, [-0.003, 0.004, 0.001])],
            &times,
        )
        .unwrap();
        for k in 0..times.len() {
            let want = [
                only1[k][0] + only2[k][0],
                only1[k][1] + only2[k][1],
                only1[k][2] + only2[k][2],
                only1[k][3] + only2[k][3],
                only1[k][4] + only2[k][4],
                only1[k][5] + only2[k][5],
            ];
            assert_close(both[k], want, &format!("sample {k}"));
        }
    }

    #[wasm_bindgen_test]
    fn coasts_hold_da_and_dix_and_rotate_de() {
        // STM structure: rows 0 (δa) and 4 (δi_x) are invariant during coasts;
        // the δe vector rotates with the secular perigee drift on an eccentric
        // chief. 200 × 30 s = 6000 s of coast after a burn at t = 0.
        let times = grid(0.0, 201, 30.0);
        let ms = vec![burn(0.0, [0.01, 0.02, -0.005])];
        let (track, _) = controlled_roe_track(&chief(), 0.0, &ms, &times).unwrap();
        let (first, last) = (track[0], track[200]);
        assert!(
            (last[0] - first[0]).abs() <= 1e-12 * first[0].abs().max(1.0),
            "δa drifted"
        );
        assert!(
            (last[4] - first[4]).abs() <= 1e-12 * first[4].abs().max(1.0),
            "δi_x drifted"
        );
        let dangle = last[3].atan2(last[2]) - first[3].atan2(first[2]);
        assert!(
            dangle.abs() > 1e-9,
            "δe must rotate under ω̇, got Δangle = {dangle}"
        );
    }

    #[wasm_bindgen_test]
    fn nonzero_t_i_epoch_matches_the_zero_epoch_case() {
        // Absolute times, durations from t_i: shifting (t_i, burn, grid) by
        // 50 000 s must reproduce the t_i = 0 result exactly (same chief epoch).
        let t_i = 50_000.0;
        let times_s: Vec<f64> = (0..11).map(|k| t_i + k as f64 * 30.0).collect();
        let (track_s, jumps_s) =
            controlled_roe_track(&chief(), t_i, &[burn(t_i, [0.01, 0.0, 0.0])], &times_s).unwrap();
        let times_0 = grid(0.0, 11, 30.0);
        let (track_0, jumps_0) =
            controlled_roe_track(&chief(), 0.0, &[burn(0.0, [0.01, 0.0, 0.0])], &times_0).unwrap();
        assert_close(jumps_s[0], jumps_0[0], "jump");
        for k in 0..11 {
            assert_close(track_s[k], track_0[k], &format!("sample {k}"));
        }
    }

    #[wasm_bindgen_test]
    fn no_maneuvers_gives_zero_track_and_no_jumps() {
        let times = grid(0.0, 11, 30.0);
        let (track, jumps) = controlled_roe_track(&chief(), 0.0, &[], &times).unwrap();
        assert!(jumps.is_empty());
        assert_eq!(track.len(), 11);
        assert!(track.iter().all(|s| *s == [0.0; 6]));
    }
}
