//! # koenig-damico-planner
//!
//! Faithful Rust re-implementation of Koenig & D'Amico's fuel-optimal impulsive
//! control algorithm (IEEE TAC 2020; see the References below).
//!
//! ## References
//!
//! Functions throughout the crate carry a `Ref:` comment citing the equation,
//! table, algorithm, or figure they implement, using these short keys:
//!
//! - **\[KD20\]** A. W. Koenig and S. D'Amico, "Fast Algorithm for Fuel-Optimal
//!   Impulsive Control of Linear Systems with Time-Varying Cost," *IEEE
//!   Transactions on Automatic Control*, 2020. DOI: 10.1109/TAC.2020.3027804
//!   (arXiv:1804.06099).
//! - **\[KGD17\]** A. W. Koenig, T. Guffanti, and S. D'Amico, "New State
//!   Transition Matrices for Spacecraft Relative Motion in Perturbed Orbits,"
//!   *Journal of Guidance, Control, and Dynamics*, 2017. DOI: 10.2514/1.G002409.
//! - **\[CD18\]** M. Chernick and S. D'Amico, "Closed-Form Optimal Impulsive
//!   Control of Spacecraft Formations Using Reachable Set Theory," AAS 18-308,
//!   2018.
//! - **\[H25\]** M. Hunter and S. D'Amico, "Fast Fuel-Optimal Constrained
//!   Impulsive Control with Application to Distributed Spacecraft," *Proc. IEEE
//!   Aerospace Conference*, 2025.

// Require every public fallible or panicking fn to document how it can fail.
#![warn(clippy::missing_errors_doc, clippy::missing_panics_doc)]

pub mod algorithm;
pub mod cost;
pub mod dynamics;
pub mod solver;
pub mod types;

// --- Core API: the entry points and results most users need. ---
pub use algorithm::{primer_history, solve, solve_from_initial_times, PrimerHistory};
pub use dynamics::Dynamics;
pub use types::{InvalidInputKind, Maneuver, PlannerError, Solution, SolveParams, TimeGrid};

// --- Problem-definition types. ---
pub use cost::{CostModel, SublevelSet};
pub use types::{Pseudostate, M, N};

// --- Convex-encoding internals (advanced use; the solve path wraps these).
// The fixed-direction QP encoding from the paper is available as
// `solver::extract_qp`; it is not re-exported at the crate root. ---
pub use solver::{min_fuel_socp, refine_socp, MinFuelSolution, RefineSolution};
pub use types::{ConicRows, Dual, FuelGenerator};
