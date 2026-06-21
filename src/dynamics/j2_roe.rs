//! J2-perturbed mean-ROE dynamics: ties the secular propagation, control-input
//! matrix `B(t)`, and state-transition matrix `Phi(t, t_f)` into the only thing
//! the algorithm needs, `Gamma(t) = Phi(t, t_f) B(t)`.

use super::b_matrix::control_input_matrix;
use super::orbit::AbsoluteOrbit;
use super::stm::state_transition;
use super::Dynamics;
use crate::types::{PlannerError, M, N};
use nalgebra::SMatrix;

/// J2 mean-ROE dynamics for a fixed chief orbit and control window `[t_i, t_f]`.
#[derive(Debug, Clone, Copy)]
pub struct J2Roe {
    chief_ti: AbsoluteOrbit,
    t_i: f64,
    t_f: f64,
}

impl J2Roe {
    /// Build from the chief's mean absolute orbit at `t_i` and the window
    /// endpoints `[t_i, t_f]` `[s]`.
    ///
    /// Ref: \[KD20\] eq. 50 / Appendix — the secular rates carry `n = sqrt(mu/a^3)`
    /// and an `a^{7/2}` denominator (real-and-finite only for a bounded ellipse:
    /// `a > 0` finite, `e ∈ [0,1)`); \[CD18\] eq. 38 / \[KD20\] B(t), p. 13 — the
    /// `1/tan(i)` term in `B(t)` is what makes the near-equatorial chief guard
    /// necessary.
    ///
    /// # Errors
    /// Returns [`PlannerError::InvalidInput`] if the chief semimajor axis is not
    /// finite and positive (`a > 0`), the chief is non-elliptic (`e ∉ [0,1)`),
    /// equatorial (`i` within `1e-9` rad of `0` or `π`, where the `tan i` term in
    /// `B(t)` is singular), or the window is not `t_f > t_i` (finite).
    pub fn new(chief_ti: AbsoluteOrbit, t_i: f64, t_f: f64) -> Result<Self, PlannerError> {
        // Bounded-ellipse preconditions on the chief (a, e together). `a` is the
        // most fundamental: `mean_motion`/`secular_rates` are NaN without a finite,
        // positive semimajor axis (\[KD20\] eq. 50), so guard it first — otherwise
        // the NaN propagates silently through `Gamma` and is mis-reported as a
        // solver failure instead of the caller-fixable input error it is.
        if !chief_ti.a.is_finite() || chief_ti.a <= 0.0 {
            return Err(PlannerError::InvalidInput(format!(
                "J2Roe: chief semimajor axis must be finite and positive (a > 0), \
                 got a = {}",
                chief_ti.a
            )));
        }
        if !(0.0..1.0).contains(&chief_ti.e) {
            return Err(PlannerError::InvalidInput(format!(
                "J2Roe: chief must be elliptic (0 <= e < 1), got e = {}",
                chief_ti.e
            )));
        }
        if chief_ti.i.sin().abs() < 1e-9 {
            return Err(PlannerError::InvalidInput(format!(
                "J2Roe: chief inclination must be bounded away from 0 and pi \
                 (B(t) has a 1/tan(i) singularity), got i = {} rad",
                chief_ti.i
            )));
        }
        if !t_i.is_finite() || !t_f.is_finite() || t_f <= t_i {
            return Err(PlannerError::InvalidInput(
                "J2Roe: window must satisfy finite t_i, t_f and t_f > t_i".into(),
            ));
        }
        Ok(Self { chief_ti, t_i, t_f })
    }
}

impl Dynamics for J2Roe {
    // Ref: [H25] eq. 2 (Gamma(t) = Phi(t, t_f) B(t)); [KD20] eq. 11 (impulsive
    // state evolution); [CD18] eq. 1 (Phi(t_f, t_k) Gamma(t_k) propagation).
    fn gamma(&self, t: f64) -> Result<SMatrix<f64, N, M>, PlannerError> {
        let orb_t = self.chief_ti.propagate(t - self.t_i);
        let orb_tf = self.chief_ti.propagate(self.t_f - self.t_i);
        let b = control_input_matrix(&orb_t)?;
        let phi = state_transition(&orb_t, &orb_tf, self.t_f - t);
        Ok(phi * b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Ref: [KD20] Table III + worked-example window [0, 117990] s.
    fn worked_example() -> J2Roe {
        let chief = AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        J2Roe::new(chief, 0.0, 117_990.0).unwrap()
    }

    #[test]
    fn new_rejects_out_of_domain_chief() {
        let ok = AbsoluteOrbit::new(25_000e3, 0.7, 40.0_f64.to_radians(), 0.0, 0.0, 0.0);
        assert!(J2Roe::new(ok, 0.0, 100.0).is_ok());

        let hyperbolic = AbsoluteOrbit::new(25_000e3, 1.2, 40.0_f64.to_radians(), 0.0, 0.0, 0.0);
        assert!(J2Roe::new(hyperbolic, 0.0, 100.0).is_err());

        let equatorial = AbsoluteOrbit::new(25_000e3, 0.7, 0.0, 0.0, 0.0, 0.0);
        assert!(J2Roe::new(equatorial, 0.0, 100.0).is_err()); // tan(i)=0 singularity

        let polar_flip = AbsoluteOrbit::new(25_000e3, 0.7, std::f64::consts::PI, 0.0, 0.0, 0.0);
        assert!(J2Roe::new(polar_flip, 0.0, 100.0).is_err()); // i = pi, tan=0

        let bad_window = AbsoluteOrbit::new(25_000e3, 0.7, 40.0_f64.to_radians(), 0.0, 0.0, 0.0);
        assert!(J2Roe::new(bad_window, 100.0, 100.0).is_err()); // t_f <= t_i
    }

    // Ref: [KD20] eq. 50 / Appendix — the secular rates carry `n = sqrt(mu/a^3)`
    // and an `a^{7/2}` denominator, real-and-finite only for a bounded ellipse
    // (`a` finite, `a > 0`). A non-positive or non-finite `a` makes `mean_motion`
    // / `secular_rates` NaN, so the chief must be rejected at the gateway with a
    // caller-fixable `InvalidInput`.
    #[test]
    fn new_rejects_nonpositive_or_nonfinite_semimajor_axis() {
        let mk = |a: f64| AbsoluteOrbit::new(a, 0.7, 40.0_f64.to_radians(), 0.0, 0.0, 0.0);
        assert!(J2Roe::new(mk(25_000e3), 0.0, 100.0).is_ok()); // baseline valid `a`
        for bad in [-1.0, 0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err =
                J2Roe::new(mk(bad), 0.0, 100.0).expect_err(&format!("a = {bad} must be rejected"));
            assert!(
                matches!(err, PlannerError::InvalidInput(_)),
                "a = {bad} must be InvalidInput, got {err:?}"
            );
        }
    }

    #[test]
    fn j2roe_is_a_dynamics_trait_object() {
        let j = worked_example();
        let _d: &dyn Dynamics = &j;
    }

    // Ref: [H25] eq. 2 (Gamma = Phi*B) with Phi(t_f, t_f) = I.
    #[test]
    fn gamma_at_tf_equals_b_since_phi_is_identity() {
        // At t = t_f, Phi(t_f, t_f) = I, so Gamma(t_f) = B(t_f).
        let j = worked_example();
        let orb_tf = j.chief_ti.propagate(j.t_f - j.t_i);
        assert_relative_eq!(
            j.gamma(j.t_f).unwrap(),
            control_input_matrix(&orb_tf).unwrap(),
            epsilon = 1e-12,
            max_relative = 1e-10
        );
    }

    // Ref: [KD20] Gamma(t) = Phi(t, t_f) B(t) (STM + B(t) displays, p. 13);
    // [H25] eq. 2 (Gamma=Phi*B), eq. 76/77 (Phi), eq. 78 (B).
    #[test]
    fn gamma_entrywise_matches_oracle() {
        let g = worked_example().gamma(16_050.0).unwrap();
        let expected = SMatrix::<f64, N, M>::from_row_slice(&[
            -4.292240669143e-04,
            4.630275430939e-04,
            0.0,
            // delta-lambda (row 2): the dominant near-degenerate entries; values are
            // of order ~1e-2 (linear-in-time drift). Cross-checked against the
            // finite-difference oracles in tests/fd_stm.rs + fd_b_matrix.rs.
            1.009859094742e-02,
            -1.131471991149e-02,
            2.136815027274e-06,
            -1.570198747958e-04,
            -2.573333198136e-05,
            -1.474305345880e-06,
            8.854647323216e-05,
            -4.013671661538e-04,
            1.991842405991e-04,
            0.0,
            0.0,
            -1.312779826620e-04,
            -2.096366596708e-06,
            5.789865377896e-06,
            -2.373367425468e-04,
        ]);
        assert_relative_eq!(g, expected, epsilon = 1e-10, max_relative = 1e-9);
    }
}
