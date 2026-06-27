//! # koenig-damico-planner
//!
//! Faithful Rust re-implementation of Koenig & D'Amico's fuel-optimal impulsive
//! control algorithm (IEEE TAC 2020; see the References below).
//!
//! ## Quick start
//!
//! ```no_run
//! use koenig_damico_planner::{solve, Pseudostate, SolveParams, TimeGrid};
//! use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
//! use koenig_damico_planner::cost::Piecewise;
//! use std::f64::consts::TAU;
//!
//! let a_c = 25_000e3; // chief semimajor axis [m], the I/O scale factor
//! let chief = AbsoluteOrbit::new(
//!     a_c, 0.7, 40f64.to_radians(), 358f64.to_radians(), 0.0, 180f64.to_radians(),
//! );
//! let dynamics = J2Roe::new(chief, 0.0, 117_990.0)?;     // validates the chief
//! let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0)?;   // validates dt > 0, t_f > t_i
//! let cost = Piecewise::new(TAU / chief.mean_motion())?; // validates period > 0
//! let w = Pseudostate::from_row_slice(&[50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]) / a_c;
//!
//! let solution = solve(&dynamics, &cost, w, grid, &SolveParams::default())?;
//! println!("{} maneuvers, total dv = {:.4} mm/s",
//!     solution.maneuvers.len(), solution.total_dv * 1e3);
//! # Ok::<(), koenig_damico_planner::PlannerError>(())
//! ```
//!
//! The concrete built-ins a caller instantiates live in the submodules:
//! `dynamics::{AbsoluteOrbit, J2Roe}` (the only [`Dynamics`] implementor) and
//! `cost::{Piecewise, Norm2, FaceMax}` (the [`CostModel`] / [`SublevelSet`]
//! implementors). A runnable version of the above is `examples/mdot.rs`
//! (`cargo run --example mdot`).
//!
//! ## Features
//!
//! - **`serde`** *(off by default)* — derives `Serialize`/`Deserialize` on the
//!   public result/wire types ([`Solution`], [`Maneuver`], [`TimeGrid`],
//!   [`SolveParams`], [`PlannerError`], [`InvalidInputKind`], and
//!   `dynamics::AbsoluteOrbit`) for the JSON request/response contract. docs.rs
//!   renders the crate with this feature enabled.
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

// --- Convex-encoding internals (advanced use; the solve path wraps these). These are reached via
// their owning modules — `solver::{min_fuel_socp, refine_socp, extract_qp, MinFuelSolution,
// RefineSolution}` and `types::{ConicRows, Dual, FuelGenerator}` — and are deliberately not
// re-exported at the crate root, to keep the root surface small. `Dual` is the exception that
// surfaces in the primary API (it types the re-exported `Solution::lambda` field and the
// `primer_history` argument); name it as `types::Dual` when you need it. ---

/// Re-export of the [`nalgebra`] version this crate's public API is built on
/// (`Pseudostate`, `Dual`, `Maneuver::dv`, `Solution::lambda`, `ConicRows`,
/// `PrimerHistory`). Use `koenig_damico_planner::nalgebra` to construct those
/// types against a matching version. A `nalgebra` *major* bump is a breaking
/// change of this crate.
pub use nalgebra;
