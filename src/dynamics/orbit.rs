//! Mean absolute Keplerian orbit, its J2 secular rates (\[KD20\] eq. 50), and
//! linear secular propagation. Mean elements in, mean elements out.

use super::constants::{J2, MU, R_E};
use super::kepler::mean_to_true;
use crate::types::PlannerError;

/// A mean absolute Keplerian orbit `[a, e, i, Omega, omega, M]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AbsoluteOrbit {
    /// Semimajor axis `[m]`.
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination `[rad]`.
    pub i: f64,
    /// Right ascension of the ascending node, Omega `[rad]`.
    pub raan: f64,
    /// Argument of perigee, omega `[rad]`.
    pub argp: f64,
    /// Mean anomaly, M `[rad]`.
    pub mean_anom: f64,
}

/// Secular rates of the slowly-varying angles under J2 (\[KD20\] eq. 50). `a`,
/// `e`, `i` are secularly constant and so have no rate.
#[derive(Debug, Clone, Copy)]
pub struct SecularRates {
    /// dOmega/dt [rad/s].
    pub raan_dot: f64,
    /// domega/dt [rad/s].
    pub argp_dot: f64,
    /// dM/dt [rad/s] (Keplerian mean motion plus the J2 secular term).
    pub mean_anom_dot: f64,
}

impl AbsoluteOrbit {
    /// Construct from the six mean elements (angles in radians).
    ///
    /// Ref: \[KD20\] mean absolute element vector `oe = [a, e, i, Omega, omega, M]`
    /// (p. 12, defined above eq. 50).
    pub fn new(a: f64, e: f64, i: f64, raan: f64, argp: f64, mean_anom: f64) -> Self {
        Self {
            a,
            e,
            i,
            raan,
            argp,
            mean_anom,
        }
    }

    /// Keplerian mean motion `n = sqrt(mu / a^3)` [rad/s].
    ///
    /// Ref: \[KGD17\] eq. 9; \[KD20\] eq. 50 (the leading Keplerian term of `Mdot`).
    pub fn mean_motion(&self) -> f64 {
        (MU / self.a.powi(3)).sqrt()
    }

    /// `eta = sqrt(1 - e^2)`.
    ///
    /// Ref: \[KGD17\] eq. 14 (the `eta` substitution).
    pub fn eta(&self) -> f64 {
        (1.0 - self.e * self.e).sqrt()
    }

    /// True anomaly `nu` `[rad]` from the current mean anomaly.
    ///
    /// # Errors
    /// Propagates [`mean_to_true`]'s errors (non-elliptic `e`).
    pub fn true_anomaly(&self) -> Result<f64, PlannerError> {
        mean_to_true(self.mean_anom, self.e)
    }

    /// J2 secular rates.
    ///
    /// Ref: \[KD20\] eq. 50 (p. 12); \[KGD17\] eq. 13 (Brouwer J2 secular rates).
    pub fn secular_rates(&self) -> SecularRates {
        let n = self.mean_motion();
        let eta = self.eta();
        let ci = self.i.cos();
        let pref = 3.0 * J2 * R_E * R_E * MU.sqrt() / self.a.powf(3.5);
        SecularRates {
            raan_dot: -pref / (2.0 * eta.powi(4)) * ci,
            argp_dot: pref / (4.0 * eta.powi(4)) * (5.0 * ci * ci - 1.0),
            mean_anom_dot: n + pref / (4.0 * eta.powi(3)) * (3.0 * ci * ci - 1.0),
        }
    }

    /// Propagate `dt` seconds: `a, e, i` constant; `Omega, omega, M` advance at
    /// their secular rates. `oe(t) = oe(t_i) + (t - t_i) * oe_dot`.
    ///
    /// Ref: \[KD20\] eq. 50 integrated (a, e, i secularly constant; angles linear
    /// in `dt`); \[KGD17\] eq. A1.
    pub fn propagate(&self, dt: f64) -> AbsoluteOrbit {
        let r = self.secular_rates();
        AbsoluteOrbit {
            a: self.a,
            e: self.e,
            i: self.i,
            raan: self.raan + r.raan_dot * dt,
            argp: self.argp + r.argp_dot * dt,
            mean_anom: self.mean_anom + r.mean_anom_dot * dt,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    // Ref: [KD20] Table III (worked-example chief mean absolute orbit).
    fn worked_example_chief() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        )
    }

    // Ref: [KD20] Table III (chief a, e -> n, eta).
    #[test]
    fn mean_motion_and_eta_match_anchors() {
        let o = worked_example_chief();
        assert_abs_diff_eq!(o.mean_motion(), 1.5971975457e-04, epsilon = 1e-13);
        assert_abs_diff_eq!(o.eta(), 0.7141428429, epsilon = 1e-9);
    }

    // Ref: [KD20] eq. 50 (evaluated at the Table III chief).
    #[test]
    fn secular_rates_match_anchors() {
        let r = worked_example_chief().secular_rates();
        assert_abs_diff_eq!(r.raan_dot, -4.9691233881e-08, epsilon = 1e-17);
        assert_abs_diff_eq!(r.argp_dot, 6.2730584504e-08, epsilon = 1e-17);
        assert_abs_diff_eq!(r.mean_anom_dot, 1.5973736883e-04, epsilon = 1e-13);
    }

    // Ref: [KD20] eq. 50 (a, e, i secularly constant; angles linear in dt).
    #[test]
    fn propagation_is_linear_and_fixes_a_e_i() {
        let o = worked_example_chief();
        let r = o.secular_rates();
        let p = o.propagate(1000.0);
        assert_abs_diff_eq!(p.a, o.a, epsilon = 1e-6);
        assert_abs_diff_eq!(p.e, o.e, epsilon = 1e-15);
        assert_abs_diff_eq!(p.i, o.i, epsilon = 1e-15);
        assert_abs_diff_eq!(p.argp, o.argp + r.argp_dot * 1000.0, epsilon = 1e-15);
        assert_abs_diff_eq!(
            p.mean_anom,
            o.mean_anom + r.mean_anom_dot * 1000.0,
            epsilon = 1e-12
        );
        let p0 = o.propagate(0.0);
        assert_abs_diff_eq!(p0.mean_anom, o.mean_anom, epsilon = 1e-15);
    }
}
