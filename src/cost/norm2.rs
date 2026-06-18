//! L2 cost `||u||_2`: the unit-ball sublevel set (Table II).

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// L2 cost `||u||_2`. Contact `g(y) = ||y||_2`, support `s(y) = y / ||y||_2`,
/// and one SOC row per time: `||Gamma^T(t) lambda||_2 <= 1`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Norm2;

impl SublevelSet for Norm2 {
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        y.norm()
    }

    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        let n = y.norm();
        if n > 0.0 {
            y / n
        } else {
            SVector::<f64, M>::zeros()
        }
    }

    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        // g(Gamma^T lambda) <= 1  <=>  one SOC row with G = Gamma^T, h = 1.
        ConicRows {
            linear: Vec::new(),
            soc: vec![(gamma_t.transpose(), 1.0)],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn contact_is_the_l2_norm() {
        // (3,4,12) has the clean norm 13 (avoids any constant-like literal).
        assert_relative_eq!(
            Norm2.contact(SVector::<f64, M>::new(3.0, 4.0, 12.0)),
            13.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            Norm2.contact(SVector::<f64, M>::new(1.0, 0.0, 0.0)),
            1.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn support_is_the_unit_direction() {
        let y = SVector::<f64, M>::new(-0.6, 0.2, -0.9);
        let s = Norm2.support(y);
        assert_relative_eq!(
            s,
            SVector::<f64, M>::new(
                -0.5454545454545454,
                0.18181818181818182,
                -0.8181818181818182
            ),
            epsilon = 1e-12
        );
        assert_relative_eq!(s.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn support_of_zero_is_zero() {
        assert_relative_eq!(
            Norm2.support(SVector::<f64, M>::zeros()),
            SVector::<f64, M>::zeros(),
            epsilon = 1e-15
        );
    }

    #[test]
    fn contact_support_identity_eq23() {
        // eq. 23: lambda . s(lambda) = g(lambda).
        for y in [
            SVector::<f64, M>::new(0.3, 0.4, 0.5),
            SVector::<f64, M>::new(-0.6, 0.2, -0.9),
        ] {
            assert_relative_eq!(y.dot(&Norm2.support(y)), Norm2.contact(y), epsilon = 1e-12);
        }
    }

    #[test]
    fn positive_homogeneity() {
        // g(a y) = a g(y) for a >= 0 (eq. 8 / Property 3).
        let y = SVector::<f64, M>::new(0.3, 0.4, 0.5);
        assert_relative_eq!(
            Norm2.contact(y * 2.5),
            2.5 * Norm2.contact(y),
            epsilon = 1e-12
        );
    }

    #[test]
    fn cone_row_matches_contact() {
        let gamma = SMatrix::<f64, N, M>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.0, 0.2, -0.3, 0.7, -0.4, 0.1,
            0.6,
        ]);
        let lam = SVector::<f64, N>::from_row_slice(&[0.5, -0.2, 0.8, 0.1, -0.6, 0.3]);
        let rows = Norm2.cone_constraints(&gamma);
        assert!(rows.linear.is_empty());
        assert_eq!(rows.soc.len(), 1);
        let (g, h) = &rows.soc[0];
        assert_relative_eq!(*h, 1.0, epsilon = 1e-15);
        // ||G lambda|| equals contact(Gamma^T lambda) equals the oracle value.
        assert_relative_eq!(
            (g * lam).norm(),
            Norm2.contact(gamma.transpose() * lam),
            epsilon = 1e-12
        );
        assert_relative_eq!((g * lam).norm(), 0.6428841264178172, epsilon = 1e-12);
    }
}
