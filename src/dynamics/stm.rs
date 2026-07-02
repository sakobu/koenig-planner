//! State-transition matrix `Phi(t, t_f)`: 6x6 quasi-nonsingular ROE
//! STM with Koenig's modified delta-lambda row 2. `dt = t_f - t`. The row-2
//! modification lives here only (the `/eta` on Phi_23/Phi_24 and the modified
//! Phi_21); `B` is unchanged. `Phi_24 = 7 kappa e_{y1} P dt / eta` is
//! intentionally nonzero (delta-lambda couples to delta-e_y under J2).

use super::constants::j2_secular_numerator;
use super::orbit::AbsoluteOrbit;
use crate::types::N;
use nalgebra::SMatrix;

/// `Phi(t, t_f)` with `dt = t_f - t`, for the chief `orb_t` at time `t`. Every
/// element is a function of `orb_t` alone: `a, e, i, omega(t)`, the mean motion,
/// and the secular `omega_dot`, with `omega(t_f) = omega(t) + omega_dot * dt`
/// derived from the eq. 50 secular drift. `a, e, i` are secularly constant, so a
/// single mean orbit fully determines `Phi(t, t_f)` — there is no second `t_f`
/// orbit to keep consistent.
///
/// Ref: \[KD20\] p. 13 — the assembled quasi-nonsingular J2 STM whose modified
/// delta-lambda row 2 (the `/eta` on `Phi_23`/`Phi_24` and the `Phi_21` drift) is
/// unique to \[KD20\] and matches neither \[KGD17\] eq. A6/A8 nor the other base
/// forms. The unmodified base STM (rows 1, 3-6) is \[KGD17\] eq. 25 / eq. A6,
/// \[CD18\] eq. 32, and \[H25\] eq. 76-77.
pub fn state_transition(orb_t: &AbsoluteOrbit, dt: f64) -> SMatrix<f64, N, N> {
    let a = orb_t.a;
    let e = orb_t.e;
    let i = orb_t.i;
    let eta = orb_t.eta();
    let n = orb_t.mean_motion();
    let w_dot = orb_t.secular_rates().argp_dot;

    let kappa = j2_secular_numerator() / (4.0 * a.powf(3.5) * eta.powi(4));
    let g = eta.powi(-2); // G = eta^-2
    let ci = i.cos();
    let p = 3.0 * ci * ci - 1.0; // P = 3 cos^2 i - 1
    let q = 5.0 * ci * ci - 1.0; // Q = 5 cos^2 i - 1
    let s = (2.0 * i).sin(); // S = sin 2i
    let t_sub = i.sin().powi(2); // sin^2 i

    let ex1 = e * orb_t.argp.cos(); // e cos omega(t)
    let ey1 = e * orb_t.argp.sin(); // e sin omega(t)

    // omega(t_f) = omega(t) + omega_dot·dt (eq. 50 secular drift); deriving it here
    // ties e_x2/e_y2 and the cos/sin(omega_dot·dt) rotation to one omega_dot·dt.
    let argp_tf = orb_t.argp + w_dot * dt;
    let ex2 = e * argp_tf.cos(); // e cos omega(t_f)
    let ey2 = e * argp_tf.sin(); // e sin omega(t_f)

    let cwd = (w_dot * dt).cos(); // cos(omega_dot dt)
    let swd = (w_dot * dt).sin(); // sin(omega_dot dt)

    let mut f = SMatrix::<f64, N, N>::zeros();
    f[(0, 0)] = 1.0;

    f[(1, 0)] = (-1.5 * n - 7.0 * kappa * eta * p) * dt;
    f[(1, 1)] = 1.0;
    f[(1, 2)] = 7.0 * kappa * ex1 * p * dt / eta;
    f[(1, 3)] = 7.0 * kappa * ey1 * p * dt / eta;
    f[(1, 4)] = -7.0 * kappa * eta * s * dt;

    f[(2, 0)] = 3.5 * kappa * ey2 * q * dt;
    f[(2, 2)] = cwd - 4.0 * kappa * ex1 * ey2 * g * q * dt;
    f[(2, 3)] = -swd - 4.0 * kappa * ey1 * ey2 * g * q * dt;
    f[(2, 4)] = 5.0 * kappa * ey2 * s * dt;

    f[(3, 0)] = -3.5 * kappa * ex2 * q * dt;
    f[(3, 2)] = swd + 4.0 * kappa * ex1 * ex2 * g * q * dt;
    f[(3, 3)] = cwd + 4.0 * kappa * ey1 * ex2 * g * q * dt;
    f[(3, 4)] = -5.0 * kappa * ex2 * s * dt;

    f[(4, 4)] = 1.0;

    f[(5, 0)] = 3.5 * kappa * s * dt;
    f[(5, 2)] = -4.0 * kappa * ex1 * g * s * dt;
    f[(5, 3)] = -4.0 * kappa * ey1 * g * s * dt;
    f[(5, 4)] = 2.0 * kappa * t_sub * dt;
    f[(5, 5)] = 1.0;

    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::AbsoluteOrbit;
    use approx::assert_relative_eq;

    fn fixture_t() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            25_000e3,
            0.3,
            50.0_f64.to_radians(),
            20.0_f64.to_radians(),
            40.0_f64.to_radians(),
            70.0_f64.to_radians(),
        )
    }

    // Ref: [KGD17] eq. 25 (Phi = I + A*tau, so Phi(tau=0) = I); [KD20] p. 13 display.
    #[test]
    fn phi_tends_to_identity_as_dt_zero() {
        let o = fixture_t();
        let phi = state_transition(&o, 0.0);
        assert_relative_eq!(phi, SMatrix::<f64, N, N>::identity(), epsilon = 1e-12);
    }

    // Ref: [KD20] Phi_24 = 7 kappa e_y1 P dt / eta (the modified row 2, p. 13).
    #[test]
    fn phi_2_4_is_nonzero() {
        // Documented: under J2 the delta-lambda row couples to delta-e_y.
        // With omega(t) = 40 deg, e_{y1} != 0, so Phi[(1,3)] != 0. Never zero it.
        let o = fixture_t();
        let phi = state_transition(&o, 39000.0);
        assert!(phi[(1, 3)].abs() > 1e-6);
    }

    // Ref: [KD20] p. 13 STM display (Phi_11..Phi_66, modified row 2); base STM
    // [KGD17] eq. 25 / eq. A6, [H25] eq. 76/77.
    #[test]
    fn entrywise_matches_oracle() {
        let o = fixture_t();
        let phi = state_transition(&o, 39000.0);
        let expected = SMatrix::<f64, N, N>::from_row_slice(&[
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            // delta-lambda secular drift is linear in dt (-1.5 n dt); this element
            // is intentionally first-order in dt. Verified entrywise by
            // tests/fd_stm.rs.
            -9.344241108678e+00,
            1.0,
            1.604820130377e-04,
            1.346603979505e-04,
            -2.612691836736e-03,
            0.0,
            2.859578380975e-04,
            0.0,
            9.999173773068e-01,
            -4.927268109377e-04,
            3.774394508970e-04,
            0.0,
            -3.404983440109e-04,
            0.0,
            5.217478630314e-04,
            1.000082372420e+00,
            -4.494281704248e-04,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.369422617739e-03,
            0.0,
            -3.952421676359e-04,
            -3.316475570890e-04,
            4.662898069914e-04,
            1.0,
        ]);
        assert_relative_eq!(phi, expected, epsilon = 1e-12, max_relative = 1e-9);
    }
}
