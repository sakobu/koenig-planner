//! Kepler's equation solve `M -> E -> nu`. Not present in any source PDF
//! (Koenig/Chernick/Hunter/ref [27] all defer to "Kepler's equation"); taken
//! from standard astrodynamics (Vallado) and verified by round-trip identity
//! and known `M -> nu` pairs, not by a PDF cross-check.

use std::f64::consts::PI;

/// Reduce an angle [rad] to the interval `[-pi, pi)`.
pub fn wrap_to_pi(x: f64) -> f64 {
    let two_pi = 2.0 * PI;
    (x + PI).rem_euclid(two_pi) - PI
}

/// Solve Kepler's equation `M = E - e sin E` for the eccentric anomaly `E` [rad].
///
/// Newton iteration with initial guess `E0 = M + e sin M`. Well-conditioned at
/// `e = 0.7` (`1 - e cos E >= 0.3`); converges in ~5-8 iterations.
pub fn mean_to_eccentric(m: f64, e: f64) -> f64 {
    let m = wrap_to_pi(m);
    let mut ecc = m + e * m.sin();
    for _ in 0..60 {
        let delta = (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
        ecc -= delta;
        if delta.abs() < 1e-14 {
            break;
        }
    }
    ecc
}

/// True anomaly `nu` [rad] from mean anomaly `M` [rad] at eccentricity `e`.
pub fn mean_to_true(m: f64, e: f64) -> f64 {
    let ecc = mean_to_eccentric(m, e);
    let eta = (1.0 - e * e).sqrt();
    (eta * ecc.sin()).atan2(ecc.cos() - e)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn wrap_reduces_to_pi_interval() {
        assert_abs_diff_eq!(wrap_to_pi(0.3), 0.3, epsilon = 1e-12);
        assert_abs_diff_eq!(wrap_to_pi(2.0 * PI + 0.3), 0.3, epsilon = 1e-12);
        assert_abs_diff_eq!(wrap_to_pi(-2.0 * PI + 0.3), 0.3, epsilon = 1e-12);
    }

    #[test]
    fn known_mean_to_true_pairs_at_e07() {
        let e = 0.7;
        assert_abs_diff_eq!(mean_to_true(0.0, e), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(mean_to_true(0.5, e), 1.9756130405, epsilon = 1e-9);
        assert_abs_diff_eq!(mean_to_true(1.0, e), 2.4310140013, epsilon = 1e-9);
        assert_abs_diff_eq!(mean_to_true(2.0, e), 2.8401081430, epsilon = 1e-9);
        // M = pi is apoapsis: nu = +/-pi. wrap_to_pi(pi) = -pi, so the solver
        // returns nu ~ -pi; compare the magnitude (both represent apoapsis).
        assert_abs_diff_eq!(mean_to_true(PI, e).abs(), PI, epsilon = 1e-9);
    }

    #[test]
    fn kepler_equation_residual_is_tiny_at_e07() {
        let e = 0.7;
        for k in 0..360 {
            let m = wrap_to_pi(k as f64 * PI / 180.0);
            let ecc = mean_to_eccentric(m, e);
            assert_abs_diff_eq!(ecc - e * ecc.sin(), m, epsilon = 1e-11);
        }
    }
}
