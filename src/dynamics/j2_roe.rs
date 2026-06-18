//! J2-perturbed mean-ROE dynamics: ties the secular propagation, control-input
//! matrix `B(t)`, and state-transition matrix `Phi(t, t_f)` into the only thing
//! the algorithm needs, `Gamma(t) = Phi(t, t_f) B(t)`.

use super::b_matrix::control_input_matrix;
use super::orbit::AbsoluteOrbit;
use super::stm::state_transition;
use super::Dynamics;
use crate::types::{M, N};
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
    /// endpoints `[t_i, t_f]` [s].
    pub fn new(chief_ti: AbsoluteOrbit, t_i: f64, t_f: f64) -> Self {
        Self { chief_ti, t_i, t_f }
    }
}

impl Dynamics for J2Roe {
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
        let orb_t = self.chief_ti.propagate(t - self.t_i);
        let orb_tf = self.chief_ti.propagate(self.t_f - self.t_i);
        let b = control_input_matrix(&orb_t);
        let phi = state_transition(&orb_t, &orb_tf, self.t_f - t);
        phi * b
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn worked_example() -> J2Roe {
        let chief = AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        J2Roe::new(chief, 0.0, 117_990.0)
    }

    #[test]
    fn j2roe_is_a_dynamics_trait_object() {
        let j = worked_example();
        let _d: &dyn Dynamics = &j;
    }

    #[test]
    fn gamma_at_tf_equals_b_since_phi_is_identity() {
        // At t = t_f, Phi(t_f, t_f) = I, so Gamma(t_f) = B(t_f).
        let j = worked_example();
        let orb_tf = j.chief_ti.propagate(j.t_f - j.t_i);
        assert_relative_eq!(
            j.gamma(j.t_f),
            control_input_matrix(&orb_tf),
            epsilon = 1e-12,
            max_relative = 1e-10
        );
    }

    #[test]
    fn gamma_entrywise_matches_oracle() {
        let g = worked_example().gamma(16_050.0);
        let expected = SMatrix::<f64, N, M>::from_row_slice(&[
            -4.292240669143e-04,
            4.630275430939e-04,
            0.0,
            1.068619416128e+03,
            -1.152778796710e+03,
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
