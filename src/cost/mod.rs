//! Cost abstractions: the unit sublevel set of the cost at a time, and the
//! time-varying selection of sublevel sets (eq. 49).

use crate::types::{ConicRows, FuelGenerator, M, N};
use nalgebra::{SMatrix, SVector};

mod private {
    /// Sealing marker: prevents [`SublevelSet`](super::SublevelSet) from being
    /// implemented outside this crate.
    ///
    /// Sealing is the *reversible* default: adding the seal is breaking (done
    /// once, pre-1.0), but *removing* it later — opening the trait if a downstream
    /// convex gauge is ever requested — is non-breaking. Until then the built-in
    /// gauges (`Norm2`, `FaceMax`, the closed [KD20] Table II set) are the only
    /// implementors, so the solver can rely on the gauge / support-function
    /// duality without guarding against incorrect external impls. `CostModel` and
    /// `Dynamics` are intentionally NOT sealed: custom time-varying cost selection
    /// and custom linear dynamics are supported extension points.
    pub trait Sealed {}
}

/// The unit sublevel set `U(1,t)` of the cost at a fixed time.
///
/// The sublevel set lives in the `M = 3` RTN control space (the same space as a
/// maneuver's Δv), so `contact` and `support` operate on `ℝ³` vectors.
///
/// This trait is **sealed**: the built-in gauges [`Norm2`] and [`FaceMax`] are
/// its only implementors. Custom downstream gauges are not supported today, but
/// the seal is reversible — it can be opened in a minor release if a downstream
/// gauge is requested. Compose the built-in gauges over time with the open
/// [`CostModel`] trait.
pub trait SublevelSet: private::Sealed {
    /// Contact function `g(y) = max_{z in U} y . z`.
    fn contact(&self, y: SVector<f64, M>) -> f64;
    /// Support direction `s(y) = argmax_{z in U} y . z`.
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M>;
    /// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1`.
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows;
    /// Primal fuel generator for the direct min-fuel SOCP: how a Δv
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
        // Piecewise carries fields, so construct it via `new`.
        let pw = Piecewise::new(39_338.811_433_158_5).unwrap();
        let _c: &dyn CostModel = &pw;
    }
}
