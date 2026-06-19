//! Time-varying piecewise cost (eq. 49): FaceMax in 2-hr perigee windows (T1),
//! Norm2 elsewhere (T2).

use super::{CostModel, FaceMax, Norm2, SublevelSet};

/// Piecewise eq.-49 selector. `T1 = { t : dist(t, nearest perigee center) < half_width }`
/// with `half_width = 1 hr` (eq. 49's 2-hr windows) and perigee centers at
/// `t_perigee0 + k·period`.
///
/// [`Piecewise::new`] assumes the chief is at apogee at `t = 0` (so the first
/// perigee is at `period/2`), which holds for the worked example (`M0 = 180°`).
/// Use [`Piecewise::with_perigee_epoch`] when the chief's perigee passage is at
/// some other time. The paper does not pin down whether `period` is the
/// Keplerian `2π/n` or the J2-perturbed period, so the caller passes the period
/// it wants (the worked example uses `2π/n` ≈ 10.93 hr).
#[derive(Debug, Clone, Copy)]
pub struct Piecewise {
    norm2: Norm2,
    facemax: FaceMax,
    period: f64,
    t_perigee0: f64,
    half_width: f64,
}

impl Piecewise {
    /// Build the eq.-49 selector for an orbit `period` [s], assuming the chief is
    /// at apogee at `t = 0` (first perigee at `period/2`). Equivalent to
    /// `with_perigee_epoch(period, period / 2.0)`. The perigee window half-width
    /// is `1 hr = 3600 s`.
    pub fn new(period: f64) -> Self {
        Self::with_perigee_epoch(period, period / 2.0)
    }

    /// Build the eq.-49 selector with an explicit perigee-passage epoch
    /// `t_perigee0` [s]; window centers are `t_perigee0 + k·period`.
    pub fn with_perigee_epoch(period: f64, t_perigee0: f64) -> Self {
        debug_assert!(
            period.is_finite() && period > 0.0,
            "Piecewise period must be finite and > 0, got {period}"
        );
        Self {
            norm2: Norm2,
            facemax: FaceMax,
            period,
            t_perigee0,
            half_width: 3600.0,
        }
    }

    /// `true` iff `t` lies within `half_width` of a perigee center
    /// `t_perigee0 + k·period`, i.e. `t` is in the eq.-49 set `T1`.
    pub fn in_perigee_window(&self, t: f64) -> bool {
        let phase = ((t - self.t_perigee0) / self.period).rem_euclid(1.0);
        let dist_frac = phase.min(1.0 - phase);
        dist_frac * self.period < self.half_width
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

    #[test]
    fn explicit_epoch_shifts_windows_off_the_apogee_default() {
        // With perigee at t=0 (not the apogee-at-0 default), the perigee passages
        // are at 0, period, 2*period, ... so t=0 is INSIDE T1 and t=period/2
        // (apogee) is OUTSIDE — the mirror image of `new`.
        let pw = Piecewise::with_perigee_epoch(40000.0, 0.0);
        assert!(pw.in_perigee_window(0.0));
        assert!(pw.in_perigee_window(40000.0));
        assert!(!pw.in_perigee_window(20000.0));
        // `new(period)` is exactly `with_perigee_epoch(period, period/2)`:
        let default = Piecewise::new(40000.0);
        assert!(!default.in_perigee_window(0.0));
        assert!(default.in_perigee_window(20000.0));
    }
}
