//! Public-API integration tests for the Phase 3 solver wrappers.

use approx::assert_relative_eq;
use koenig_planner::cost::{FaceMax, Norm2, SublevelSet};
use koenig_planner::{extract_qp, refine_socp};
use nalgebra::{SMatrix, SVector};

const N: usize = 6;
const M: usize = 3;

fn gamma_top_identity() -> SMatrix<f64, N, M> {
    let mut g = SMatrix::<f64, N, M>::zeros();
    for i in 0..M {
        g[(i, i)] = 1.0;
    }
    g
}
fn gamma_bottom_identity() -> SMatrix<f64, N, M> {
    let mut g = SMatrix::<f64, N, M>::zeros();
    for i in 0..M {
        g[(M + i, i)] = 1.0;
    }
    g
}

#[test]
fn refine_then_extract_hands_off_through_public_api() {
    // Mixed problem: FaceMax on l1..3 + Norm2 (SOC) on l4..6, w=(0,0,1,0,0,1).
    // c* = sqrt(3)+1 (validated in unit tests); use it as the QP budget.
    let rows = vec![
        FaceMax.cone_constraints(&gamma_top_identity()),
        Norm2.cone_constraints(&gamma_bottom_identity()),
    ];
    let w = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]);

    let refined = refine_socp(&w, &rows).unwrap();
    assert_relative_eq!(refined.objective, 3.0_f64.sqrt() + 1.0, epsilon = 1e-6);
    assert!(refined.objective >= 0.0);

    // Hand off to the QP: two directions that exactly reconstruct w
    // (y1 = e3, y2 = e6); budget = c* is slack, so alpha = (1, 1), residual 0.
    let y1 = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
    let y2 = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);
    let q = SMatrix::<f64, N, N>::identity();
    let alpha = extract_qp(&w, &[y1, y2], &q, refined.objective).unwrap();
    assert_relative_eq!(alpha[0], 1.0, epsilon = 1e-6);
    assert_relative_eq!(alpha[1], 1.0, epsilon = 1e-6);

    let werr = w - alpha[0] * y1 - alpha[1] * y2;
    assert!(werr.norm() < 1e-6);
}
