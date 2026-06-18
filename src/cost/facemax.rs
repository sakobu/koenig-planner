//! Face-max cost `max(V_face u)` for the tetrahedral fixed-attitude occulter
//! (eq. 47-48). Implemented in Phase 2.

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// Face-max cost `max_k y^T w_k` over `W = [0, V_vertex]`. Linear rows per time.
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceMax;

impl SublevelSet for FaceMax {
    #[allow(unused_variables)]
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        unimplemented!("Phase 2: g(y) = max_k y^T w_k")
    }
    #[allow(unused_variables)]
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        unimplemented!("Phase 2: s(y) = argmax_k column")
    }
    #[allow(unused_variables)]
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        unimplemented!("Phase 2/3: linear rows w_k^T Gamma^T lambda <= 1 for all k")
    }
}
