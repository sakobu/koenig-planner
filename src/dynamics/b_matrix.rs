//! Control-input matrix `B(t)`: GVE map from an RTN Delta-v [m/s] to
//! a dimensionless mean-ROE change. Columns = R, T, N thrust on the deputy.
//! `theta = omega + nu`; `nu` is the true anomaly from `M` via Kepler.

use super::constants::MU;
use super::orbit::AbsoluteOrbit;
use crate::types::{PlannerError, M, N};
use nalgebra::SMatrix;

/// `B(t)` evaluated at `orbit`, including the `sqrt(a/mu)` scaling. The `[B_ij]`
/// block depends only on `e, i, omega, nu`; `a` enters solely through the scale.
///
/// # Errors
/// Propagates [`AbsoluteOrbit::true_anomaly`]'s errors (non-elliptic `e`).
pub fn control_input_matrix(orbit: &AbsoluteOrbit) -> Result<SMatrix<f64, N, M>, PlannerError> {
    let e = orbit.e;
    let i = orbit.i;
    let argp = orbit.argp;
    let eta = orbit.eta();
    let nu = orbit.true_anomaly()?;
    let theta = argp + nu;
    let ecn = 1.0 + e * nu.cos(); // 1 + e cos nu
    let tan_i = i.tan();
    let scale = (orbit.a / MU).sqrt();

    let mut b = SMatrix::<f64, N, M>::zeros();
    b[(0, 0)] = (2.0 / eta) * e * nu.sin();
    b[(0, 1)] = (2.0 / eta) * ecn;
    b[(1, 0)] = -2.0 * eta * eta / ecn;
    b[(2, 0)] = eta * theta.sin();
    b[(2, 1)] = eta * ((2.0 + e * nu.cos()) * theta.cos() + e * argp.cos()) / ecn;
    b[(2, 2)] = eta * e * argp.sin() * theta.sin() / (tan_i * ecn);
    b[(3, 0)] = -eta * theta.cos();
    b[(3, 1)] = eta * ((2.0 + e * nu.cos()) * theta.sin() + e * argp.sin()) / ecn;
    b[(3, 2)] = -eta * e * argp.cos() * theta.sin() / (tan_i * ecn);
    b[(4, 2)] = eta * theta.cos() / ecn;
    b[(5, 2)] = eta * theta.sin() / ecn;
    Ok(b * scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::AbsoluteOrbit;
    use approx::assert_relative_eq;

    fn fixture() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            25_000e3,
            0.3,
            50.0_f64.to_radians(),
            20.0_f64.to_radians(),
            40.0_f64.to_radians(),
            70.0_f64.to_radians(),
        )
    }

    #[test]
    fn zero_structure_is_correct() {
        let b = control_input_matrix(&fixture()).unwrap();
        assert_eq!(b[(0, 2)], 0.0);
        assert_eq!(b[(1, 1)], 0.0);
        assert_eq!(b[(1, 2)], 0.0);
        assert_eq!(b[(4, 0)], 0.0);
        assert_eq!(b[(4, 1)], 0.0);
        assert_eq!(b[(5, 0)], 0.0);
        assert_eq!(b[(5, 1)], 0.0);
    }

    #[test]
    fn entrywise_matches_oracle() {
        let b = control_input_matrix(&fixture()).unwrap();
        let expected = SMatrix::<f64, N, M>::from_row_slice(&[
            1.523378438764e-04,
            4.849959064280e-04,
            0.0,
            -4.934524718319e-04,
            0.0,
            0.0,
            1.379310012159e-04,
            -3.468028551194e-04,
            2.416221553306e-05,
            1.950635808382e-04,
            3.371317264229e-04,
            -2.879540716656e-05,
            0.0,
            0.0,
            -2.111780503734e-04,
            0.0,
            0.0,
            1.493256701105e-04,
        ]);
        assert_relative_eq!(b, expected, epsilon = 1e-12, max_relative = 1e-9);
    }

    #[test]
    fn b_scales_as_sqrt_a_over_mu() {
        // [B_ij] is a-independent; only the sqrt(a/mu) scale carries a.
        // Quadrupling a (same e,i,omega,M) scales B by sqrt(4) = 2.
        let o1 = fixture();
        let o4 = AbsoluteOrbit::new(4.0 * o1.a, o1.e, o1.i, o1.raan, o1.argp, o1.mean_anom);
        let b1 = control_input_matrix(&o1).unwrap();
        let b4 = control_input_matrix(&o4).unwrap();
        assert_relative_eq!(b4, b1 * 2.0, epsilon = 1e-12, max_relative = 1e-10);
    }
}
