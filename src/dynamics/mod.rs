//! Dynamics abstraction. The algorithm only ever needs `Gamma(t) = Phi(t,t_f) B(t)`.

use crate::types::{PlannerError, M, N};
use nalgebra::SMatrix;

/// Maps an impulse at time `t` into pseudostate space via `Gamma(t) = Phi(t,t_f) B(t)`.
pub trait Dynamics {
    /// `Gamma(t)` in R^{6x3}: pseudostate change per unit Delta-v `[m/s]` applied at `t` `[s]`.
    ///
    /// # Errors
    /// Returns [`PlannerError`] if evaluating `B(t)` requires an out-of-domain
    /// Kepler solve (non-elliptic chief). For a `J2Roe` built via its validating
    /// constructor this cannot occur, but the trait is fallible so other
    /// implementations may report domain failures rather than panic.
    fn gamma(&self, t: f64) -> Result<SMatrix<f64, N, M>, PlannerError>;
}

pub mod b_matrix;
pub mod constants;
pub mod j2_roe;
pub mod kepler;
pub mod orbit;
pub mod stm;
pub use j2_roe::J2Roe;
pub use orbit::AbsoluteOrbit;
