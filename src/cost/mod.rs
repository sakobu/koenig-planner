//! Cost abstractions: the unit sublevel set of the cost at a time, and the
//! time-varying selection of sublevel sets (eq. 49).

use crate::types::{ConicRows, FuelGenerator, M, N};
use nalgebra::{SMatrix, SVector};

/// The unit sublevel set `U(1,t)` of the cost at a fixed time.
pub trait SublevelSet {
    /// Contact function `g(y) = max_{z in U} y . z`.
    fn contact(&self, y: SVector<f64, M>) -> f64;
    /// Support direction `s(y) = argmax_{z in U} y . z`.
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M>;
    /// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1`.
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows;
    /// Primal fuel generator for the direct min-fuel SOCP (Phase 5b): how a Δv
    /// in this sublevel set is built from solver variables and charged.
    fn fuel_generator(&self) -> FuelGenerator;
}

/// Time-varying cost = piecewise selection of sublevel sets (eq. 49).
pub trait CostModel {
    /// The sublevel set active at time `t`.
    fn at(&self, t: f64) -> &dyn SublevelSet;
}

pub mod facemax;
pub mod norm2;
pub mod piecewise;

pub use facemax::FaceMax;
pub use norm2::Norm2;
pub use piecewise::Piecewise;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cost_types_wire_to_their_traits() {
        let _s: &dyn SublevelSet = &Norm2;
        let _f: &dyn SublevelSet = &FaceMax;
        // Piecewise now carries fields, so construct it via `new`.
        let pw = Piecewise::new(39_338.811_433_158_5);
        let _c: &dyn CostModel = &pw;
    }
}
