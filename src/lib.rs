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

pub use algorithm::solve;
pub use cost::{CostModel, SublevelSet};
pub use dynamics::Dynamics;
pub use solver::{extract_qp, refine_socp, RefineSolution};
pub use types::{
    ConicRows, Dual, FuelGenerator, Maneuver, PlannerError, Pseudostate, Solution, SolveParams,
    TimeGrid, M, N,
};
