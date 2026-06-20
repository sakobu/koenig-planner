//! Face-max cost `max(V_face u)` for the tetrahedral fixed-attitude occulter
//! (eq. 47-48). The algorithm consumes the cost only through its contact /
//! support / cone-constraint forms, all expressed via `W = [0 | V_vertex]`
//! (Table II); `V_face` itself is the cost definition and is not needed here.

use std::sync::LazyLock;

use super::SublevelSet;
use crate::types::{ConicRows, FuelGenerator, M, N};
use nalgebra::{SMatrix, SVector};

/// Face-max cost. Contact `g(y) = max(0, max_k y . v_k)` over the four
/// `V_vertex` columns (the origin column of `W` supplies the `max(0, .)`);
/// support `s(y)` is the argmax column (origin -> zero when all `y . v_k <= 0`);
/// linear rows `(Gamma v_k)^T lambda <= 1` for each vertex column.
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceMax;

/// The four tetrahedral `V_vertex` support directions (eq. 48), computed once.
static VERTEX_COLUMNS: LazyLock<[SVector<f64, M>; 4]> = LazyLock::new(|| {
    let a = (2.0_f64 / 3.0).sqrt();
    let b = (1.0_f64 / 3.0).sqrt();
    [
        SVector::<f64, M>::new(a, 0.0, -b),
        SVector::<f64, M>::new(-a, 0.0, -b),
        SVector::<f64, M>::new(0.0, a, b),
        SVector::<f64, M>::new(0.0, -a, b),
    ]
});

fn vertex_columns() -> [SVector<f64, M>; 4] {
    *VERTEX_COLUMNS
}

impl SublevelSet for FaceMax {
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        // g(y) = max over W = [0, V_vertex]; the origin column contributes 0.0,
        // so g(y) >= 0 always (and g(0) = 0).
        vertex_columns()
            .iter()
            .map(|v| y.dot(v))
            .fold(0.0, f64::max)
    }

    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        // argmax column of W; ties resolve to the lowest index; the origin
        // column (zero vector, value 0) wins when every y . v_k <= 0.
        let mut best = SVector::<f64, M>::zeros();
        let mut best_val = 0.0_f64;
        for v in vertex_columns() {
            let val = y.dot(&v);
            if val > best_val {
                best_val = val;
                best = v;
            }
        }
        best
    }

    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        // g(Gamma^T lambda) <= 1  <=>  (Gamma v_k)^T lambda <= 1 for each vertex
        // column (the origin column gives the vacuous 0 <= 1 and is omitted).
        let linear = vertex_columns()
            .iter()
            .map(|v| (gamma_t * v, 1.0))
            .collect();
        ConicRows {
            linear,
            soc: Vec::new(),
        }
    }

    fn fuel_generator(&self) -> FuelGenerator {
        // Unit ball U(1) = conv{0, V_vertex columns}; its gauge is
        //   f(v) = min{ Σₖ θₖ : Σₖ θₖ vₖ = v, θ ≥ 0 }.
        // (The eq. 48 V_face matrix is the cost-defining form and is not used here;
        // the algorithm's geometry is fully carried by the V_vertex columns above.)
        FuelGenerator::Polytope(vertex_columns().to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn vertices_are_unit_and_tetrahedral() {
        let cols = vertex_columns();
        for v in &cols {
            assert_relative_eq!(v.norm(), 1.0, epsilon = 1e-12);
        }
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    assert_relative_eq!(cols[i].dot(&cols[j]), -1.0 / 3.0, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn contact_known_directions() {
        let s23 = (2.0_f64 / 3.0).sqrt();
        let s13 = (1.0_f64 / 3.0).sqrt();
        assert_relative_eq!(
            FaceMax.contact(SVector::<f64, M>::new(1.0, 0.0, 0.0)),
            s23,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            FaceMax.contact(SVector::<f64, M>::new(0.0, 0.0, 1.0)),
            s13,
            epsilon = 1e-12
        );
    }

    #[test]
    fn contact_of_zero_is_zero() {
        assert_relative_eq!(
            FaceMax.contact(SVector::<f64, M>::zeros()),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn support_is_argmax_vertex() {
        let cols = vertex_columns();
        // y = (1,0,0): unique argmax = vertex column 0.
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::new(1.0, 0.0, 0.0)),
            cols[0],
            epsilon = 1e-12
        );
        // y = (0,0,1): tie between columns 2 and 3; lowest index (2) wins.
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::new(0.0, 0.0, 1.0)),
            cols[2],
            epsilon = 1e-12
        );
        // y = (-0.6,0.2,-0.9): argmax = vertex column 1.
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::new(-0.6, 0.2, -0.9)),
            cols[1],
            epsilon = 1e-12
        );
    }

    #[test]
    fn support_of_zero_is_zero() {
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::zeros()),
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
            assert_relative_eq!(
                y.dot(&FaceMax.support(y)),
                FaceMax.contact(y),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn positive_homogeneity() {
        // g(a y) = a g(y) for a >= 0 (eq. 8 / Property 3).
        let y = SVector::<f64, M>::new(0.3, 0.4, 0.5);
        assert_relative_eq!(
            FaceMax.contact(y * 2.5),
            2.5 * FaceMax.contact(y),
            epsilon = 1e-12
        );
    }

    #[test]
    fn vertex_face_transcription_cross_check() {
        // V_face (eq. 48), 4x3. Used only to cross-check the transcription of
        // both matrices: f(v_k) = max_row(V_face v_k) is the same (= 1/9) for
        // every vertex column v_k of V_vertex.
        let s23 = (2.0_f64 / 3.0).sqrt();
        let s13 = (1.0_f64 / 3.0).sqrt();
        let v_face = SMatrix::<f64, 4, 3>::from_row_slice(&[
            -s23, 0.0, s13, s23, 0.0, s13, 0.0, -s23, -s13, 0.0, s23, -s13,
        ]) / 3.0;
        let cols = vertex_columns();
        for v in &cols {
            assert_relative_eq!((v_face * v).max(), 1.0 / 9.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn fuel_generator_is_polytope_of_unit_vertices() {
        use crate::types::FuelGenerator;
        match FaceMax.fuel_generator() {
            FuelGenerator::Polytope(dirs) => {
                assert_eq!(dirs.len(), 4);
                // Same four unit tetrahedral directions used by contact/support.
                for (d, v) in dirs.iter().zip(vertex_columns()) {
                    assert_relative_eq!(*d, v, epsilon = 1e-12);
                    assert_relative_eq!(d.norm(), 1.0, epsilon = 1e-12);
                }
            }
            other => panic!("expected Polytope, got {other:?}"),
        }
    }

    #[test]
    fn cone_rows_match_contact() {
        let gamma = SMatrix::<f64, N, M>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.0, 0.2, -0.3, 0.7, -0.4, 0.1,
            0.6,
        ]);
        let lam = SVector::<f64, N>::from_row_slice(&[0.5, -0.2, 0.8, 0.1, -0.6, 0.3]);
        let rows = FaceMax.cone_constraints(&gamma);
        assert!(rows.soc.is_empty());
        assert_eq!(rows.linear.len(), 4);
        for (_, b) in &rows.linear {
            assert_relative_eq!(*b, 1.0, epsilon = 1e-15);
        }
        let max_row = rows
            .linear
            .iter()
            .map(|(a, _)| a.dot(&lam))
            .fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(
            max_row,
            FaceMax.contact(gamma.transpose() * lam),
            epsilon = 1e-12
        );
        assert_relative_eq!(max_row, 0.37230594560185404, epsilon = 1e-12);
    }
}
