//! Time-varying piecewise cost (eq. 49): FaceMax in 2-hr perigee windows (T1),
//! Norm2 elsewhere (T2).

use super::{CostModel, FaceMax, Norm2, SublevelSet};

/// Piecewise eq.-49 selector. `T1 = { t : |t - (k+0.5) period| < half_width }`
/// with `half_width = 1 hr` (eq. 49's 2-hr windows). The centers
/// `(k+0.5) period` land on perigee for the worked example because its chief
/// starts at apogee (`M0 = 180 deg`) at `t = 0`. The paper does not pin down
/// whether `period` is the Keplerian `2 pi / n` or the J2-perturbed period, so
/// the caller passes the period it wants (the worked example uses `2 pi / n`,
/// approx 10.93 hr, consistent with the paper's rounded 10.92 hr).
#[derive(Debug, Clone, Copy)]
pub struct Piecewise {
    norm2: Norm2,
    facemax: FaceMax,
    period: f64,
    half_width: f64,
}

impl Piecewise {
    /// Build the eq.-49 selector for an orbit period `period` [s]; the perigee
    /// window half-width is `1 hr = 3600 s`.
    pub fn new(period: f64) -> Self {
        Self {
            norm2: Norm2,
            facemax: FaceMax,
            period,
            half_width: 3600.0,
        }
    }

    /// `true` iff `t` lies within `half_width` of a perigee center
    /// `(k+0.5) period`, i.e. `t` is in the eq.-49 set `T1`.
    pub fn in_perigee_window(&self, t: f64) -> bool {
        let frac = (t / self.period).rem_euclid(1.0);
        (frac - 0.5).abs() * self.period < self.half_width
    }
}

impl CostModel for Piecewise {
    fn at(&self, t: f64) -> &dyn SublevelSet {
        if self.in_perigee_window(t) {
            &self.facemax
        } else {
            &self.norm2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::M;
    use approx::assert_relative_eq;
    use nalgebra::SVector;

    #[test]
    fn perigee_window_boundaries() {
        // period 40000 s: first perigee center at 20000, half_width 3600, so
        // T1 is the open interval (16400, 23600). Probe 1 s inside / outside
        // each edge. The exact 3600 s boundary is excluded by the strict `<`,
        // but it is a floating-point knife-edge, so it is not asserted directly.
        let pw = Piecewise::new(40000.0);
        // Center, and 3599 s either side -> inside T1.
        assert!(pw.in_perigee_window(20000.0));
        assert!(pw.in_perigee_window(16401.0));
        assert!(pw.in_perigee_window(23599.0));
        // 3601 s either side -> outside T1.
        assert!(!pw.in_perigee_window(16399.0));
        assert!(!pw.in_perigee_window(23601.0));
        // Apogees (t = 0 and t = period) -> outside T1.
        assert!(!pw.in_perigee_window(0.0));
        assert!(!pw.in_perigee_window(40000.0));
        // Second orbit's perigee at (1.5) * 40000 = 60000 -> inside T1.
        assert!(pw.in_perigee_window(60000.0));
    }

    #[test]
    fn at_selects_facemax_in_window_norm2_outside() {
        let pw = Piecewise::new(40000.0);
        let ex = SVector::<f64, M>::new(1.0, 0.0, 0.0);
        // Inside T1 -> FaceMax: g(ex) = sqrt(2/3).
        assert_relative_eq!(
            pw.at(20000.0).contact(ex),
            (2.0_f64 / 3.0).sqrt(),
            epsilon = 1e-12
        );
        // Outside T1 -> Norm2: g(ex) = 1.
        assert_relative_eq!(pw.at(0.0).contact(ex), 1.0, epsilon = 1e-12);
    }
}
