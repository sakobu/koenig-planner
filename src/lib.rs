//! # koenig-planner
//!
//! Faithful Rust reimplementation of Koenig & D'Amico's fuel-optimal impulsive
//! control algorithm (IEEE TAC 2020). See `docs/Planner.pdf` and
//! `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`.

pub mod cost;
pub mod dynamics;
pub mod types;

pub use cost::{CostModel, SublevelSet};
pub use dynamics::Dynamics;
pub use types::{
    ConicRows, Dual, Maneuver, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N,
};
