//! # koenig-planner
//!
//! Faithful Rust re-implementation of Koenig & D'Amico's fuel-optimal impulsive
//! control algorithm (IEEE TAC 2020). See `docs/Planner.pdf`.

pub mod algorithm;
pub mod cost;
pub mod dynamics;
pub mod solver;
pub mod types;

// --- Core API: the entry points and results most users need. ---
pub use algorithm::{solve, solve_from_initial_times};
pub use dynamics::Dynamics;
pub use types::{Maneuver, PlannerError, Solution, SolveParams, TimeGrid};

// --- Problem-definition types. ---
pub use cost::{CostModel, SublevelSet};
pub use types::{Pseudostate, M, N};

// --- Convex-encoding internals (advanced use; the solve path wraps these).
// The fixed-direction QP encoding from the paper is available as
// `solver::extract_qp`; it is not re-exported at the crate root. ---
pub use solver::{min_fuel_socp, refine_socp, MinFuelSolution, RefineSolution};
pub use types::{ConicRows, Dual, FuelGenerator};
