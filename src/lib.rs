//! # koenig-planner
//!
//! Faithful Rust re-implementation of Koenig & D'Amico's fuel-optimal impulsive
//! control algorithm (IEEE TAC 2020). See `docs/Planner.pdf` and
//! `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`.

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
// `extract_qp` is intentionally NOT re-exported here: it is the paper's
// superseded fixed-direction QP, retained under `solver::extract_qp`. ---
pub use solver::{min_fuel_socp, refine_socp, MinFuelSolution, RefineSolution};
pub use types::{ConicRows, Dual, FuelGenerator};
