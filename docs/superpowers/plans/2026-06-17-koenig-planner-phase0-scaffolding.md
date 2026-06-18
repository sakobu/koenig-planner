# Koenig Planner — Phase 0 (Scaffolding) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the `koenig-planner` Rust library crate — module tree, core value types, trait seams, stub impls, binaries/examples, and CI — so the whole skeleton compiles and `cargo test` is green on real (non-stub) tests, ready for Phase 1 to fill in the dynamics.

**Architecture:** A single library crate (Approach A from the design §4) rooted at the **project root** (`Cargo.toml` and `src/` sit next to the existing `docs/` so the spec's `docs/Planner.pdf` references stay valid). Trait seams (`Dynamics`, `SublevelSet`, `CostModel`) mirror the paper's math abstractions; concrete impls (`J2Roe`, `Norm2`, `FaceMax`, `Piecewise`) and the solver/algorithm layers are present as compiling stubs that `unimplemented!()` with a phase pointer. The only behavior implemented in Phase 0 is the value types (`TimeGrid`, `SolveParams`, `Maneuver`, …), which carry the phase's real tests.

**Tech Stack:** Rust 2021 (toolchain ≥ 1.92); `nalgebra` 0.35 (static-dim linear algebra); `clarabel` 0.11 (conic solver, declared now, used in Phase 3); `thiserror` 2.0 (errors); `approx` 0.5 (dev, float-tolerant asserts); `csv` 1.4 + `plotters` 0.3 (optional, behind a `validation` feature, used in Phase 6). GitHub Actions for CI.

**Verification status:** This entire scaffold was built in a scratch crate and run through the exact CI gates before this plan was written — `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test` (default **and** `--all-features`) all pass; 10 tests green; bins/examples compile and run; `clarabel`/`nalgebra` build on Rust 1.92. Every code block below is transcribed from that verified build (already rustfmt-ordered).

## Global Constraints

These apply to every task; copy values verbatim.

- **Crate location:** repo root. `Cargo.toml` at `/Users/sarkismelkonian/Desktop/Koenig Planner/Cargo.toml`; sources under `src/`. Do **not** nest under a `koenig-planner/` subdirectory — the spec references `docs/Planner.pdf` relative to the crate root.
- **Package:** `name = "koenig-planner"`, `version = "0.0.0"`, `edition = "2021"`, `rust-version = "1.92"`, `license = "MIT OR Apache-2.0"`.
- **Dimensions:** `pub const N: usize = 6;` (ROE state), `pub const M: usize = 3;` (RTN Δv). Never hardcode `6`/`3` where `N`/`M` read clearer.
- **Dependency pins (exact):** `nalgebra = "0.35"`, `clarabel = "0.11"`, `thiserror = "2.0"`; dev `approx = "0.5"`; optional `csv = "1.4"` + `plotters = "0.3"` behind `validation = ["dep:csv", "dep:plotters"]`.
- **Stubs:** non-Phase-0 trait methods/functions bodies are `unimplemented!("Phase N: …")` with a phase pointer; annotate the fn with `#[allow(unused_variables)]` so unused params don't trip `-D warnings`. Stub types `#[derive(Debug, Clone, Copy, Default)]` (unit structs).
- **Lint gate (must stay green):** `cargo clippy --all-targets --all-features -- -D warnings`. Practical consequences: prefix unused bindings with `_`; provide `is_empty` whenever a public `len` exists; compare floats in tests with `approx`, never `assert_eq!` on `f64`.
- **Formatting:** rustfmt default. Imports are sorted case-sensitively — note `Solution` sorts **before** `SolveParams` (`l` < `v`). Run `cargo fmt --all` before every commit.
- **Frame/units (forward context, no Phase-0 code):** RTN (≡ RIC) frame; mean-elements-in/mean-elements-out; `a_c` applied once at the I/O boundary, never baked into `B`/`Φ`. Documented in §5.4–5.5 of the spec for Phase 1.
- **Commit discipline:** one commit per task, conventional-commit message. `Cargo.lock` **is** committed (the crate ships binaries). `/target` and `.DS_Store` are git-ignored.

---

## File Structure

| File | Responsibility | Introduced in |
|---|---|---|
| `Cargo.toml` | package metadata, dependency manifest, `validation` feature | Task 1 |
| `.gitignore` | ignore `/target`, `.DS_Store` | Task 1 |
| `src/lib.rs` | crate docs, module declarations, public re-exports | grows Tasks 1–4 |
| `src/types.rs` | `N`, `M`, `Pseudostate`, `Dual`, `Maneuver`, `TimeGrid`, `SolveParams`, `Solution`, `ConicRows`, `PlannerError` + unit tests | Task 1 |
| `src/dynamics/mod.rs` | `Dynamics` trait | Task 2 |
| `src/dynamics/j2_roe.rs` | `J2Roe` stub (Phase 1 target) | Task 2 |
| `src/cost/mod.rs` | `SublevelSet` + `CostModel` traits | Task 3 |
| `src/cost/norm2.rs` | `Norm2` stub | Task 3 |
| `src/cost/facemax.rs` | `FaceMax` stub | Task 3 |
| `src/cost/piecewise.rs` | `Piecewise` stub (CostModel) | Task 3 |
| `src/solver/mod.rs` | solver module root | Task 4 |
| `src/solver/refine_socp.rs` | refinement SOCP (Phase 3 target; doc-only) | Task 4 |
| `src/solver/extract_qp.rs` | extraction QP (Phase 3 target; doc-only) | Task 4 |
| `src/algorithm/mod.rs` | `solve(...)` public entry-point signature | Task 4 |
| `src/algorithm/{init,refine,extract}.rs` | Algorithms 1/2/3 (Phase 4 targets; doc-only) | Task 4 |
| `tests/api.rs` | public-API smoke tests | Task 4 |
| `examples/mdot.rs` | worked-example binary (Phase 5 target; stub) | Task 5 |
| `src/bin/monte_carlo.rs` | Monte Carlo harness (Phase 6 target; stub) | Task 5 |
| `.github/workflows/ci.yml` | fmt + clippy + build + test gate | Task 5 |

---

## Task 1: Bootstrap crate & core value types

Establishes the git repo + crate, locks the dependency manifest, and implements `src/types.rs` — the only file with real Phase-0 behavior, so it carries the phase's first green tests (grid counts, default params, error display, conic-rows default).

**Files:**
- Create: `Cargo.toml`
- Create: `.gitignore`
- Create: `src/lib.rs`
- Create: `src/types.rs` (incl. `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: nothing (entry task).
- Produces (relied on by every later task):
  - `pub const N: usize = 6; pub const M: usize = 3;`
  - `pub type Pseudostate = SVector<f64, N>;` and `pub type Dual = SVector<f64, N>;`
  - `pub struct Maneuver { pub t: f64, pub dv: SVector<f64, M> }`
  - `pub struct TimeGrid { pub t_i: f64, pub t_f: f64, pub dt: f64 }` with `uniform(t_i, t_f, dt) -> Self`, `len(&self) -> usize`, `is_empty(&self) -> bool`, `time(&self, idx: usize) -> f64`, `times(&self) -> impl Iterator<Item = f64> + '_`
  - `pub struct SolveParams { pub n_coarse: usize, pub n_init: usize, pub eps_cost: f64, pub eps_remove: f64, pub q: SMatrix<f64, N, N> }` with `Default` = Table III
  - `pub struct Solution { pub maneuvers: Vec<Maneuver>, pub total_dv: f64, pub iterations: usize, pub residual: f64, pub lambda: Dual }`
  - `pub struct ConicRows { pub linear: Vec<(SVector<f64, N>, f64)>, pub soc: Vec<(SMatrix<f64, M, N>, f64)> }` (`Default`)
  - `pub enum PlannerError { SolverFailed(String), NotConverged{max_iters,achieved,target}, KeplerDivergence{m,e}, InvalidInput(String) }`

- [ ] **Step 1: Initialize the git repo and crate skeleton**

```bash
cd "/Users/sarkismelkonian/Desktop/Koenig Planner"
git init
cargo init --lib --vcs none .
rm -f src/lib.rs   # replaced below
```

- [ ] **Step 2: Write `Cargo.toml`**

```toml
[package]
name = "koenig-planner"
version = "0.0.0"
edition = "2021"
rust-version = "1.92"
license = "MIT OR Apache-2.0"
description = "Faithful Rust reimplementation of the Koenig-D'Amico fuel-optimal impulsive control algorithm (IEEE TAC 2020)"

[dependencies]
nalgebra = "0.35"
clarabel = "0.11"
thiserror = "2.0"

[dependencies.csv]
version = "1.4"
optional = true

[dependencies.plotters]
version = "0.3"
optional = true

[dev-dependencies]
approx = "0.5"

[features]
validation = ["dep:csv", "dep:plotters"]
```

- [ ] **Step 3: Write `.gitignore`**

```gitignore
/target
.DS_Store
```

- [ ] **Step 4: Write `src/lib.rs` (types-only at this task; later tasks add modules)**

```rust
//! # koenig-planner
//!
//! Faithful Rust reimplementation of Koenig & D'Amico's fuel-optimal impulsive
//! control algorithm (IEEE TAC 2020). See `docs/Planner.pdf` and
//! `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`.

pub mod types;

pub use types::{
    ConicRows, Dual, Maneuver, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N,
};
```

- [ ] **Step 5: Write the failing tests inside `src/types.rs`**

Add this module at the bottom of `src/types.rs` (the `use super::*` items resolve once Step 7 lands):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn dimensions_match_spec() {
        assert_eq!(N, 6);
        assert_eq!(M, 3);
    }

    #[test]
    fn worked_example_grid_has_3934_times() {
        let g = TimeGrid::uniform(0.0, 117990.0, 30.0);
        assert_eq!(g.len(), 3934);
        assert_abs_diff_eq!(g.time(0), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(g.time(g.len() - 1), 117990.0, epsilon = 1e-6);
        assert_eq!(g.times().count(), 3934);
    }

    #[test]
    fn hunter_grid_has_3901_times() {
        let g = TimeGrid::uniform(0.0, 39000.0, 10.0);
        assert_eq!(g.len(), 3901);
    }

    #[test]
    fn default_params_match_table_iii() {
        let p = SolveParams::default();
        assert_eq!(p.n_coarse, 20);
        assert_eq!(p.n_init, 6);
        assert_abs_diff_eq!(p.eps_cost, 0.01, epsilon = 1e-12);
        assert_abs_diff_eq!(p.eps_remove, 0.01, epsilon = 1e-12);
        assert_eq!(p.q, SMatrix::<f64, N, N>::identity());
    }

    #[test]
    fn conic_rows_default_is_empty() {
        let c = ConicRows::default();
        assert!(c.linear.is_empty());
        assert!(c.soc.is_empty());
    }
}
```

- [ ] **Step 6: Run the tests to verify they fail (red)**

Run: `cargo test --lib`
Expected: FAIL to compile — `cannot find type TimeGrid` / `SolveParams` / `ConicRows` in this scope (the types in Step 7 don't exist yet).

- [ ] **Step 7: Implement the types above the test module in `src/types.rs`**

Put this at the **top** of `src/types.rs`, before the `#[cfg(test)]` block from Step 5:

```rust
//! Core value types: dimensions, pseudostate/maneuver, the uniform time grid,
//! solver parameters, the solver result, the conic-row placeholder, and errors.

use nalgebra::{SMatrix, SVector};
use thiserror::Error;

/// State dimension: 6 quasi-nonsingular relative orbital elements (ROEs).
pub const N: usize = 6;

/// Control dimension: 3 RTN Delta-v components (R, T, N).
pub const M: usize = 3;

/// A pseudostate / ROE vector in R^6 (dimensionless unless scaled by `a_c`).
pub type Pseudostate = SVector<f64, N>;

/// The dual variable lambda in R^6 (outward reachable-set normal).
pub type Dual = SVector<f64, N>;

/// An impulsive maneuver: a Delta-v [m/s] in the RTN frame applied at time `t` [s].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Maneuver {
    /// Application time [s], measured from `t_i`.
    pub t: f64,
    /// Delta-v [m/s], RTN components (R, T, N).
    pub dv: SVector<f64, M>,
}

/// A uniform, endpoint-inclusive time grid over `[t_i, t_f]` with step `dt`.
///
/// The worked example (Table III) uses a 30 s grid over `[0, 117990]` -> 3934
/// candidate times; the Hunter cross-check uses a 10 s grid over `[0, 39000]`
/// -> 3901 candidate times.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeGrid {
    /// Initial time `t_i` [s].
    pub t_i: f64,
    /// Final time `t_f` [s].
    pub t_f: f64,
    /// Grid step `dt` [s].
    pub dt: f64,
}

impl TimeGrid {
    /// Build a uniform grid with step `dt` over `[t_i, t_f]`.
    pub fn uniform(t_i: f64, t_f: f64, dt: f64) -> Self {
        Self { t_i, t_f, dt }
    }

    /// Number of grid points, inclusive of both endpoints.
    pub fn len(&self) -> usize {
        ((self.t_f - self.t_i) / self.dt).round() as usize + 1
    }

    /// A grid always has at least one point; provided for lint-completeness.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The time [s] at grid index `idx`.
    pub fn time(&self, idx: usize) -> f64 {
        self.t_i + (idx as f64) * self.dt
    }

    /// Iterator over all grid times [s].
    pub fn times(&self) -> impl Iterator<Item = f64> + '_ {
        (0..self.len()).map(move |i| self.time(i))
    }
}

/// Tunable parameters for the three-step algorithm (Table III defaults).
#[derive(Debug, Clone)]
pub struct SolveParams {
    /// Coarse-sample count `|T^d|` for Algorithm 1 initialization (Table III: 20).
    pub n_coarse: usize,
    /// Initial candidate-time count `n_init` (Table III: 6).
    pub n_init: usize,
    /// Convergence tolerance `eps_cost` (Table III: 0.01).
    pub eps_cost: f64,
    /// Slack-removal tolerance `eps_remove` (Table III: 0.01).
    pub eps_remove: f64,
    /// Positive-definite weight `Q` for the extraction QP (identity in the example).
    pub q: SMatrix<f64, N, N>,
}

impl Default for SolveParams {
    fn default() -> Self {
        Self {
            n_coarse: 20,
            n_init: 6,
            eps_cost: 0.01,
            eps_remove: 0.01,
            q: SMatrix::<f64, N, N>::identity(),
        }
    }
}

/// The planner output: the maneuver set plus convergence diagnostics.
#[derive(Debug, Clone)]
pub struct Solution {
    /// Maneuvers `{t, Delta-v}`, one per optimal time in `T^opt`.
    pub maneuvers: Vec<Maneuver>,
    /// Total fuel cost sum of ||Delta-v_j|| [m/s].
    pub total_dv: f64,
    /// Algorithm 2 iteration count.
    pub iterations: usize,
    /// Relative residual ||w_err|| / ||w||.
    pub residual: f64,
    /// Optimal dual lambda_opt.
    pub lambda: Dual,
}

/// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1` for one candidate time.
///
/// The concrete representation is finalized in Phase 3 (solver layer); this
/// placeholder exists so the [`crate::SublevelSet`] trait and its impls compile
/// in Phase 0. Linear rows encode `a^T lambda <= b`; SOC rows encode
/// `||G lambda||_2 <= h`.
#[derive(Debug, Clone, Default)]
pub struct ConicRows {
    /// Linear rows `(a, b)` with `a^T lambda <= b`.
    pub linear: Vec<(SVector<f64, N>, f64)>,
    /// Second-order-cone rows `(G, h)` with `||G lambda||_2 <= h`.
    pub soc: Vec<(SMatrix<f64, M, N>, f64)>,
}

/// Errors surfaced by the planner.
#[derive(Debug, Error)]
pub enum PlannerError {
    /// The convex (SOCP/QP) solver failed.
    #[error("convex solver failed: {0}")]
    SolverFailed(String),
    /// Algorithm 2 did not reach `max_t g <= 1 + eps_cost` within the iteration cap.
    #[error("refinement did not converge in {max_iters} iterations (max_t g = {achieved}, target <= {target})")]
    NotConverged {
        /// Iteration cap that was hit.
        max_iters: usize,
        /// Achieved `max_t g`.
        achieved: f64,
        /// Target threshold `1 + eps_cost`.
        target: f64,
    },
    /// The Kepler Newton iteration failed to converge.
    #[error("Kepler solve diverged for M = {m} rad, e = {e}")]
    KeplerDivergence {
        /// Mean anomaly [rad].
        m: f64,
        /// Eccentricity.
        e: f64,
    },
    /// An input precondition was violated.
    #[error("invalid input: {0}")]
    InvalidInput(String),
}
```

> **Note on `is_empty`:** clippy's `len_without_is_empty` requires it whenever a public `len` exists, hence the always-`false` `is_empty`. The `as usize` cast is intentional and clippy-clean under default lints.

- [ ] **Step 8: Run the tests to verify they pass (green)**

Run: `cargo test --lib`
Expected: PASS — `5 passed; 0 failed` (`dimensions_match_spec`, `worked_example_grid_has_3934_times`, `hunter_grid_has_3901_times`, `default_params_match_table_iii`, `conic_rows_default_is_empty`).

- [ ] **Step 9: Format, lint, then commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat: bootstrap koenig-planner crate and core value types"
```
Expected: clippy clean; commit created.

---

## Task 2: Dynamics trait + `J2Roe` stub

Adds the `dynamics` module: the `Dynamics` trait (the only thing Algorithms 1–3 need from the model) and the `J2Roe` stub whose `gamma` lands in Phase 1.

**Files:**
- Create: `src/dynamics/mod.rs`
- Create: `src/dynamics/j2_roe.rs`
- Modify: `src/lib.rs` (add `pub mod dynamics;` + re-export)
- Test: `src/dynamics/j2_roe.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `types::{M, N}`.
- Produces:
  - `pub trait Dynamics { fn gamma(&self, t: f64) -> SMatrix<f64, N, M>; }`
  - `pub struct J2Roe;` implementing `Dynamics` (stub) — re-exported as `dynamics::J2Roe` and `crate::Dynamics`.

- [ ] **Step 1: Write the failing test in `src/dynamics/j2_roe.rs`**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j2roe_is_a_dynamics_trait_object() {
        // Phase 0 wiring check: J2Roe constructs and is object-safe as Dynamics.
        let _d: &dyn Dynamics = &J2Roe;
    }
}
```

- [ ] **Step 2: Run it to verify it fails (red)**

Run: `cargo test --lib`
Expected: FAIL to compile — `cannot find type J2Roe` / `cannot find trait Dynamics`.

- [ ] **Step 3: Write `src/dynamics/mod.rs`**

```rust
//! Dynamics abstraction. The algorithm only ever needs `Gamma(t) = Phi(t,t_f) B(t)`.

use crate::types::{M, N};
use nalgebra::SMatrix;

/// Maps an impulse at time `t` into pseudostate space via `Gamma(t) = Phi(t,t_f) B(t)`.
pub trait Dynamics {
    /// `Gamma(t)` in R^{6x3}: pseudostate change per unit Delta-v [m/s] applied at `t` [s].
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M>;
}

pub mod j2_roe;
pub use j2_roe::J2Roe;
```

- [ ] **Step 4: Write `src/dynamics/j2_roe.rs` (above the test module from Step 1)**

```rust
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
```

- [ ] **Step 5: Wire the module into `src/lib.rs`**

Change the module/re-export block so it reads (note alphabetical module order and the `Dynamics` re-export):

```rust
pub mod dynamics;
pub mod types;

pub use dynamics::Dynamics;
pub use types::{
    ConicRows, Dual, Maneuver, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N,
};
```

- [ ] **Step 6: Run the test to verify it passes (green)**

Run: `cargo test --lib`
Expected: PASS — now `6 passed; 0 failed` (5 from Task 1 + `j2roe_is_a_dynamics_trait_object`).

- [ ] **Step 7: Format, lint, then commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat: add Dynamics trait and J2Roe stub"
```
Expected: clippy clean; commit created.

---

## Task 3: Cost traits + `Norm2`/`FaceMax`/`Piecewise` stubs

Adds the `cost` module: the `SublevelSet` and `CostModel` traits plus the three stub cost types Phase 2 will implement.

**Files:**
- Create: `src/cost/mod.rs`
- Create: `src/cost/norm2.rs`
- Create: `src/cost/facemax.rs`
- Create: `src/cost/piecewise.rs`
- Modify: `src/lib.rs` (add `pub mod cost;` + re-exports)
- Test: `src/cost/mod.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `types::{ConicRows, M, N}`.
- Produces:
  - `pub trait SublevelSet { fn contact(&self, y: SVector<f64,M>) -> f64; fn support(&self, y: SVector<f64,M>) -> SVector<f64,M>; fn cone_constraints(&self, gamma_t: &SMatrix<f64,N,M>) -> ConicRows; }`
  - `pub trait CostModel { fn at(&self, t: f64) -> &dyn SublevelSet; }`
  - `pub struct Norm2;` + `pub struct FaceMax;` implementing `SublevelSet` (stubs)
  - `pub struct Piecewise;` implementing `CostModel` (stub)
  - re-exported as `crate::{SublevelSet, CostModel}` and `cost::{Norm2, FaceMax, Piecewise}`.

- [ ] **Step 1: Write the failing test in `src/cost/mod.rs`**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cost_types_wire_to_their_traits() {
        let _s: &dyn SublevelSet = &Norm2;
        let _f: &dyn SublevelSet = &FaceMax;
        let _c: &dyn CostModel = &Piecewise;
    }
}
```

- [ ] **Step 2: Run it to verify it fails (red)**

Run: `cargo test --lib`
Expected: FAIL to compile — `cannot find type Norm2` / `FaceMax` / `Piecewise` / traits not in scope.

- [ ] **Step 3: Write `src/cost/mod.rs` (above the test module from Step 1)**

```rust
//! Cost abstractions: the unit sublevel set of the cost at a time, and the
//! time-varying selection of sublevel sets (eq. 49).

use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// The unit sublevel set `U(1,t)` of the cost at a fixed time.
pub trait SublevelSet {
    /// Contact function `g(y) = max_{z in U} y . z`.
    fn contact(&self, y: SVector<f64, M>) -> f64;
    /// Support direction `s(y) = argmax_{z in U} y . z`.
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M>;
    /// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1`.
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows;
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
```

- [ ] **Step 4: Write `src/cost/norm2.rs`**

```rust
//! L2 cost `||u||_2`: the unit-ball sublevel set. Implemented in Phase 2.

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// L2 cost `||u||_2`. One SOC row per time: `||Gamma^T(t) lambda||_2 <= 1`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Norm2;

impl SublevelSet for Norm2 {
    #[allow(unused_variables)]
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        unimplemented!("Phase 2: g(y) = ||y||_2")
    }
    #[allow(unused_variables)]
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        unimplemented!("Phase 2: s(y) = y / ||y||_2")
    }
    #[allow(unused_variables)]
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        unimplemented!("Phase 2/3: one SOC row ||Gamma^T lambda||_2 <= 1")
    }
}
```

- [ ] **Step 5: Write `src/cost/facemax.rs`**

```rust
//! Face-max cost `max(V_face u)` for the tetrahedral fixed-attitude occulter
//! (eq. 47-48). Implemented in Phase 2.

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// Face-max cost `max_k y^T w_k` over `W = [0, V_vertex]`. Linear rows per time.
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceMax;

impl SublevelSet for FaceMax {
    #[allow(unused_variables)]
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        unimplemented!("Phase 2: g(y) = max_k y^T w_k")
    }
    #[allow(unused_variables)]
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        unimplemented!("Phase 2: s(y) = argmax_k column")
    }
    #[allow(unused_variables)]
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        unimplemented!("Phase 2/3: linear rows w_k^T Gamma^T lambda <= 1 for all k")
    }
}
```

- [ ] **Step 6: Write `src/cost/piecewise.rs`**

```rust
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
```

- [ ] **Step 7: Wire the module into `src/lib.rs`**

```rust
pub mod cost;
pub mod dynamics;
pub mod types;

pub use cost::{CostModel, SublevelSet};
pub use dynamics::Dynamics;
pub use types::{
    ConicRows, Dual, Maneuver, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N,
};
```

- [ ] **Step 8: Run the test to verify it passes (green)**

Run: `cargo test --lib`
Expected: PASS — `7 passed; 0 failed`.

- [ ] **Step 9: Format, lint, then commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat: add SublevelSet/CostModel traits and Norm2/FaceMax/Piecewise stubs"
```
Expected: clippy clean; commit created.

---

## Task 4: Solver + algorithm skeleton, public `solve`, integration tests

Completes the library surface: the `solver` and `algorithm` modules (doc-only stubs for the Phase 3/4 internals), the public `solve(...)` entry-point signature that pins what Phases 4–5 must produce, and the `tests/api.rs` integration smoke tests.

**Files:**
- Create: `src/solver/mod.rs`
- Create: `src/solver/refine_socp.rs`
- Create: `src/solver/extract_qp.rs`
- Create: `src/algorithm/mod.rs`
- Create: `src/algorithm/init.rs`
- Create: `src/algorithm/refine.rs`
- Create: `src/algorithm/extract.rs`
- Create: `tests/api.rs`
- Modify: `src/lib.rs` (add `pub mod solver; pub mod algorithm;` + `pub use algorithm::solve;`)

**Interfaces:**
- Consumes: `cost::CostModel`, `dynamics::Dynamics`, `types::{PlannerError, Pseudostate, Solution, SolveParams, TimeGrid}`.
- Produces:
  - `pub fn solve<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C, w: Pseudostate, grid: TimeGrid, params: &SolveParams) -> Result<Solution, PlannerError>` (stub) — re-exported as `crate::solve`.
  - `solver::{refine_socp, extract_qp}` and `algorithm::{init, refine, extract}` module placeholders.

- [ ] **Step 1: Write the failing integration test `tests/api.rs`**

```rust
//! Public-API smoke tests for the Phase 0 scaffold.

use approx::assert_abs_diff_eq;
use koenig_planner::{Maneuver, PlannerError, SolveParams, TimeGrid, M, N};
use nalgebra::SVector;

#[test]
fn reexports_are_reachable() {
    assert_eq!(N, 6);
    assert_eq!(M, 3);
}

#[test]
fn maneuver_constructs_and_exposes_fields() {
    let m = Maneuver {
        t: 16050.0,
        dv: SVector::<f64, 3>::new(9.68e-3, -23.02e-3, -25.56e-3),
    };
    assert_abs_diff_eq!(m.t, 16050.0, epsilon = 1e-9);
    assert_eq!(m.dv.len(), 3);
}

#[test]
fn default_params_are_table_iii() {
    let p = SolveParams::default();
    assert_eq!(p.n_init, 6);
    assert_eq!(p.n_coarse, 20);
}

#[test]
fn error_displays_message() {
    let e = PlannerError::InvalidInput("bad w".into());
    assert!(e.to_string().contains("bad w"));
}

#[test]
fn worked_and_hunter_grid_counts() {
    assert_eq!(TimeGrid::uniform(0.0, 117990.0, 30.0).len(), 3934);
    assert_eq!(TimeGrid::uniform(0.0, 39000.0, 10.0).len(), 3901);
}
```

- [ ] **Step 2: Run it to verify it fails (red)**

Run: `cargo test --test api`
Expected: FAIL — at this point the integration test compiles against existing types and **passes**, *except* it cannot exercise `solve` yet; the real red gate is the next step's `cargo build` of `lib.rs` once it references the not-yet-created `algorithm` module. (If you prefer a hard red here, temporarily add `let _ = koenig_planner::solve::<koenig_planner::dynamics::J2Roe, koenig_planner::cost::Piecewise>;` — it fails to resolve until Step 5, then remove it.)

> Practical note: these five assertions are real and independent of `solve`; they stay green. The purpose of Task 4 is to add the `solve` signature + module skeleton **without breaking** them, so treat "all green after wiring" as the gate.

- [ ] **Step 3: Write the `solver` module (root + doc-only stubs)**

`src/solver/mod.rs`:

```rust
//! Convex-solver wrappers around `clarabel`: the refinement SOCP (eq. 40) and
//! the extraction QP (Algorithm 3). Implemented in Phase 3.

pub mod extract_qp;
pub mod refine_socp;
```

`src/solver/refine_socp.rs`:

```rust
//! Builds and solves eq. 40 over a candidate-time set (linear + SOC cones from
//! each time's `cone_constraints`), maximize -> minimize. Implemented in Phase 3.
```

`src/solver/extract_qp.rs`:

```rust
//! The Algorithm 3 QP: solve for nonnegative magnitudes that minimize the
//! weighted pseudostate residual. Implemented in Phase 3.
```

- [ ] **Step 4: Write the `algorithm` module (root with `solve` + doc-only stubs)**

`src/algorithm/mod.rs`:

```rust
//! Orchestration of the three-step algorithm (Init -> Refine -> Extract).

mod extract;
mod init;
mod refine;

use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{PlannerError, Pseudostate, Solution, SolveParams, TimeGrid};

/// Solve the fuel-optimal impulsive control problem.
///
/// Wires Algorithm 1 (init) -> Algorithm 2 (refine) -> Algorithm 3 (extract)
/// with `Gamma(t)` caching. Implemented in Phases 4-5.
#[allow(unused_variables)]
pub fn solve<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
) -> Result<Solution, PlannerError> {
    unimplemented!("Phases 4-5 wire init -> refine -> extract")
}
```

`src/algorithm/init.rs`:

```rust
//! Algorithm 1 - Initialization: pick the `n_init` coarse times with the
//! largest contact value as `T^est`. Implemented in Phase 4.
```

`src/algorithm/refine.rs`:

```rust
//! Algorithm 2 - Iterative Refinement: solve eq. 40 on `T^est`, drop slack
//! times, add violated local maxima, until convergence. Implemented in Phase 4.
```

`src/algorithm/extract.rs`:

```rust
//! Algorithm 3 - Control-Input Extraction: optimal directions, then the QP for
//! magnitudes. Implemented in Phase 4.
```

- [ ] **Step 5: Finalize `src/lib.rs` (full module tree + re-exports)**

```rust
//! # koenig-planner
//!
//! Faithful Rust reimplementation of Koenig & D'Amico's fuel-optimal impulsive
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
pub use types::{
    ConicRows, Dual, Maneuver, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N,
};
```

- [ ] **Step 6: Run the full test suite to verify green**

Run: `cargo test`
Expected: PASS — lib unit tests `7 passed` + `tests/api.rs` `5 passed`; bins/examples none yet.

- [ ] **Step 7: Format, lint, then commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat: add solver/algorithm skeleton, public solve signature, API smoke tests"
```
Expected: clippy clean; commit created.

---

## Task 5: Binaries/examples skeleton + CI + full-gate verification

Adds the remaining build targets (the worked-example binary and the Monte Carlo harness, as stubs), the GitHub Actions CI workflow, and runs the **exact** CI gates locally — the Phase 0 exit criterion.

**Files:**
- Create: `examples/mdot.rs`
- Create: `src/bin/monte_carlo.rs`
- Create: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: `crate::{N, M}` (the example prints them to prove linkage).
- Produces: two runnable targets (`cargo run --example mdot`, `cargo run --bin monte_carlo`) and a CI pipeline running fmt/clippy/build/test.

- [ ] **Step 1: Write `examples/mdot.rs`**

```rust
//! Worked example (Table III -> Table IV, Fig. 7 data). Implemented in Phase 5.

fn main() {
    println!(
        "koenig-planner scaffold - worked example pending (Phase 5). N={}, M={}",
        koenig_planner::N,
        koenig_planner::M
    );
}
```

- [ ] **Step 2: Write `src/bin/monte_carlo.rs`**

```rust
//! Monte Carlo harness (Fig. 8 / Fig. 9). Implemented in Phase 6.

fn main() {
    println!("koenig-planner Monte Carlo harness pending (Phase 6).");
}
```

- [ ] **Step 3: Verify both targets compile and run**

Run: `cargo run --quiet --example mdot && cargo run --quiet --bin monte_carlo`
Expected:
```
koenig-planner scaffold - worked example pending (Phase 5). N=6, M=3
koenig-planner Monte Carlo harness pending (Phase 6).
```

- [ ] **Step 4: Write `.github/workflows/ci.yml`**

```yaml
name: CI

on:
  push:
  pull_request:

env:
  CARGO_TERM_COLOR: always

jobs:
  check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v5
      # plotters (validation feature, pulled by --all-features) needs system
      # fontconfig on Linux; absent on the runner, so install it. macOS uses
      # CoreText and does not need this — a local gate cannot catch the gap.
      - name: Install system deps (fontconfig, for the plotters validation feature)
        run: sudo apt-get update && sudo apt-get install -y libfontconfig1-dev
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: rustfmt, clippy
      - uses: Swatinem/rust-cache@v2
      - name: Format
        run: cargo fmt --all -- --check
      - name: Clippy
        run: cargo clippy --all-targets --all-features -- -D warnings
      - name: Build
        run: cargo build --all-targets --all-features
      - name: Test
        run: cargo test --all-targets --all-features
```

- [ ] **Step 5: Run the exact CI gates locally (the Phase 0 exit criterion)**

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo build --all-targets --all-features
cargo test --all-targets --all-features
```
Expected: all four succeed; `cargo test` reports `10 passed` (5 lib `types` + 1 `dynamics` + 1 `cost` are lib unit tests = `7 passed` in the lib binary, plus `tests/api.rs` `5 passed`). Wait for the actual counts: lib unit = 7, integration = 5 → 12 test executions across binaries; bins/examples contribute 0 tests but must compile.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "ci: add binaries/examples skeleton and GitHub Actions pipeline"
```
Expected: commit created.

- [ ] **Step 7: (Optional, to literally "run CI") connect a GitHub remote**

CI executes only once the repo is pushed to GitHub. If/when desired:

```bash
gh repo create koenig-planner --private --source=. --remote=origin --push
```
Expected: the `CI` workflow appears under the repo's Actions tab and runs the four gates from Step 5. Until a remote exists, Step 5's local run is the binding gate.

---

## Self-Review

**1. Spec coverage (Phase 0 line in design §6 + architecture §4 + types in §4.1):**

| Spec requirement | Task |
|---|---|
| Cargo lib + deps (`nalgebra`, `clarabel`, `thiserror`, `approx`, `csv`/`plotters`) | Task 1 (manifest) |
| `types.rs`: `Pseudostate = SVector<f64,6>` | Task 1 |
| `Maneuver{t, dv:SVector<f64,3>}` | Task 1 |
| `TimeGrid` | Task 1 |
| `SolveParams` | Task 1 |
| error enum | Task 1 (`PlannerError`) |
| Empty trait defs compile (`Dynamics`) | Task 2 |
| Empty trait defs compile (`SublevelSet`, `CostModel`) | Task 3 |
| `ConicRows` (referenced by `cone_constraints`) | Task 1 (placeholder) |
| Full module tree (`dynamics/cost/solver/algorithm` + `j2_roe/norm2/facemax/piecewise/refine_socp/extract_qp/init/refine/extract`) | Tasks 2–4 |
| `examples/mdot.rs`, `src/bin/monte_carlo.rs` | Task 5 |
| **Exit:** `cargo test` green on stubs | Tasks 1–4 (10 tests) |
| **Exit:** CI runs | Task 5 |

No Phase-0 spec item is unaddressed. (`solve` signature and `Solution`/`Dual` types are added beyond the literal list to pin the public API early — additive, not a gap.)

**2. Placeholder scan:** No `TODO`/`fill in`/"add error handling"/"similar to Task N" left in the plan. Stub bodies are deliberate, complete `unimplemented!("Phase N: …")` calls (the phase's defined output), each with the real signature shown — not plan placeholders.

**3. Type consistency:** `N`/`M`, `Pseudostate`, `Dual`, `Maneuver`, `TimeGrid`, `SolveParams`, `Solution`, `ConicRows`, `PlannerError` are defined once in Task 1 and referenced with identical names/signatures in Tasks 2–5. `Dynamics::gamma`, `SublevelSet::{contact,support,cone_constraints}`, `CostModel::at`, and `solve<D,C>(…)` match between their definition task and the `tests/api.rs`/lib re-exports. Import orderings shown are rustfmt-stable (`Solution` before `SolveParams`; modules alphabetical).

One correction folded in during review: Task 4 Step 2's "red" is a compile-of-`lib` gate rather than a failing assertion, because the five `tests/api.rs` assertions are genuinely independent of `solve` and stay green — the note there says so explicitly so the executor isn't surprised when the test is green before `solve` exists.
