//! L2 cost `||u||_2`: the unit-ball sublevel set. Implemented in Phase 2.

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// L2 cost `||u||_2`. One SOC row per time: `||Gamma^T(t) lambda||_2 <= 1`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Norm2;

impl SublevelSet for Norm2 {
    #[allow(unused_variables)]
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        unimplemented!("Phase 2: g(y) = ||y||_2")
    }
    #[allow(unused_variables)]
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        unimplemented!("Phase 2: s(y) = y / ||y||_2")
    }
    #[allow(unused_variables)]
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        unimplemented!("Phase 2/3: one SOC row ||Gamma^T lambda||_2 <= 1")
    }
}
