//! J2-perturbed mean-ROE dynamics (Appendix). Implemented in Phase 1.

use super::Dynamics;
use crate::types::{M, N};
use nalgebra::SMatrix;

/// J2 mean-ROE dynamics: mean-element secular propagation, `B(t)`, ROE STM
/// `Phi(t,t_f)`, and `Gamma(t) = Phi B`. Fields and construction land in Phase 1.
#[derive(Debug, Clone, Copy, Default)]
pub struct J2Roe;

impl Dynamics for J2Roe {
    #[allow(unused_variables)]
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
        unimplemented!("Phase 1: J2 mean-ROE Gamma(t) = Phi(t,t_f) B(t)")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j2roe_is_a_dynamics_trait_object() {
        // Phase 0 wiring check: J2Roe constructs and is object-safe as Dynamics.
        let _d: &dyn Dynamics = &J2Roe;
    }
}
