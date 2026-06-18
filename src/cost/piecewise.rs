//! Time-varying piecewise cost (eq. 49): FaceMax in 2-hr perigee windows (T1),
//! Norm2 elsewhere (T2). Implemented in Phase 2.

use super::{CostModel, SublevelSet};

/// Piecewise eq.-49 selector. Holds the two sublevel sets and the window
/// geometry (orbit period, perigee offset); fields land in Phase 2.
#[derive(Debug, Clone, Copy, Default)]
pub struct Piecewise;

impl CostModel for Piecewise {
    #[allow(unused_variables)]
    fn at(&self, t: f64) -> &dyn SublevelSet {
        unimplemented!("Phase 2: select FaceMax (T1) vs Norm2 (T2) per eq. 49 windows")
    }
}
