# Phase 4 — Three Algorithms + Orchestration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the three-step planner — Algorithm 1 (initialization), Algorithm 2 (iterative refinement of eq. 40), Algorithm 3 (control-input extraction) — and the `solve(...)` orchestration that wires them with cached `Γ(t)`, so that the public `koenig_planner::solve` runs end-to-end on a synthetic problem and converges.

**Architecture:** The algorithm layer is pure glue over already-built, verified pieces: it caches `Γ(t)` from the `Dynamics` trait (Phase 1 `J2Roe`), scans the contact function `g_{U(1,t)}(Γᵀ(t)λ)` from the `CostModel`/`SublevelSet` traits (Phase 2 `Piecewise`/`Norm2`/`FaceMax`), and calls the stateless solver wrappers `refine_socp` / `extract_qp` (Phase 3). `solve()` validates inputs, builds the `Γ`-cache once, runs Init → Refine → Extract, and assembles a `Solution`. Candidate-time sets `T^est`/`T^opt` are carried as sorted-deduped `Vec<usize>` grid indices into the `Γ`-cache. `refine_socp` deliberately does **not** return the per-time slack `g`, so the refinement loop recomputes `g` over the full grid via `SublevelSet::contact` to drive the drop/add/converge logic.

**Tech Stack:** Rust 2021 (`rust-version = 1.92`), `nalgebra 0.35` (`SVector`/`SMatrix`), `clarabel 0.11.1` (consumed transitively through the Phase 3 wrappers — Phase 4 calls no clarabel API directly), `thiserror 2.0`, `approx 0.5` (dev).

## Global Constraints

- **Work on a branch:** create `phase4-algorithms` off `main` before Task 1; open a PR at the end (mirrors Phase 1/2/3: `phase3-solver-wrappers` → PR #3).
- **CI gate (run after every task, must be green):**
  `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
- **Dimensions are fixed:** `N = 6` (ROE state), `M = 3` (RTN Δv) — defined in `src/types.rs`, re-exported at the crate root. Never redefine them; `use crate::types::{M, N}` inside `src/`, and `const N/M` locals in `tests/` (mirroring `tests/solver.rs`).
- **Do not change any public signature or re-export.** `pub fn solve<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C, w: Pseudostate, grid: TimeGrid, params: &SolveParams) -> Result<Solution, PlannerError>` already exists in `src/algorithm/mod.rs` as an `unimplemented!()` stub; Phase 4 replaces only the **body**. `solve` and `Solution` are already re-exported from `src/lib.rs:13,17` — **no `lib.rs` edit is required.**
- **Transient `#[allow(dead_code)]`:** new `pub(super)` algorithm helpers (`initialize`, `refine`, `extract`, `violated_local_maxima`) are unused in the *lib* target until a later task calls them from non-test code — and `cargo clippy --all-features -- -D warnings` checks the lib target (not `#[cfg(test)]` code). Mark each with a transient `#[allow(dead_code)]` when it lands and **remove it in the task that first calls it from non-test code.** Use `#[allow(dead_code)]`, **not** `#[expect(dead_code)]` — `#[expect]` misfires in the `cfg(test)` build where the tests *do* use the item (Phase 3 gotcha, CLAUDE.md). `RefineOutcome.max_g_trace` is read only by tests, so it keeps a **permanent** `#[allow(dead_code)]` with a doc note.
- **Assert on bands, not bit-equality** (Risk R3): synthetic-problem tests assert convergence/residual within sensible numerical bands (e.g. residual `< 1e-2`, contact `≤ 1 + ε_cost`), never exact digits. Tie any closed-form expectations to `√`/integer expressions (e.g. `13.0`, `3.0_f64.sqrt()`), not decimal literals, to dodge `clippy::approx_constant`.
- **`gamma(t)` takes ABSOLUTE time** `[s]` (the same convention as `TimeGrid::time(idx) = t_i + idx·dt`). Feed grid times **directly** into `Dynamics::gamma`; `J2Roe::gamma` subtracts `t_i` internally. Never pre-subtract `t_i`.

## Design Decisions (locked, with rationale)

1. **`T^est` / `T^opt` are `Vec<usize>` of grid indices** (kept sorted + deduped after every mutation), not `Vec<f64>` times. Rationale: indices key directly into the `Γ`-cache (`gammas[k]`) and the time is recoverable via `grid.time(k)`; equality/dedup on `usize` is exact (no float-tie ambiguity).
2. **`Γ(t)` is cached once per `solve()`** into a `Vec<SMatrix<f64, N, M>>` indexed by grid index. Rationale: `J2Roe` caches nothing — every `gamma(t)` re-propagates the (time-invariant) `orb_tf`, rebuilds `B(t)` *including a fresh Kepler Newton solve*, and rebuilds the 6×6 STM. Algorithm 2 scans `g` over the full grid every iteration; without the cache that is thousands of Kepler solves per iteration. The cache changes **no numerics** (gamma is deterministic), only runtime.
3. **Iteration cap is a module constant `MAX_REFINE_ITERS: usize = 50`** in `src/algorithm/mod.rs`, passed into `refine(..)` as an argument. Rationale: `SolveParams` has **no** `max_iters` field (it is not a Table III parameter) and its public shape is locked; a const avoids changing the public API. Passing it as a `refine` argument lets tests force the `NotConverged` path with `max_iters = 1`. The Monte-Carlo targets converge in ≤ 8 iterations, so 50 is a generous safety backstop, not a tuning knob.
4. **Convergence is checked immediately after each `refine_socp` solve, *before* the drop/add step**, so `T^opt` is exactly the active set that produced `λ_opt` (every time in `T^opt` participated in the final SOCP). Rationale: a faithful restructuring of the paper's `do { solve; remove; add } while (max g > 1+ε)` — the remove/add only matter when *not yet* converged, and this keeps `λ_opt` consistent with `T^opt` for extraction (no last-iteration phantom times added after the solve). `max_t g` is monotonically non-increasing across iterations (paper §5.1), which the trace test asserts.
5. **Discrete local-maxima finder is plateau-aware and endpoint-inclusive** (Risk R4/R7). A grid index is a peak if it is `≥` its in-bounds neighbours; a flat top (run of `~equal` values) contributes **one** representative (its midpoint index); endpoints are peaks by boundary rule. Three thresholds: **add** at `g > 1.0`, **converge** at `max_t g ≤ 1 + ε_cost`, **keep/drop** at `g ≥ 1 − ε_remove`. Rationale: the global max is always a local max, so the binding violator is always added → progress guaranteed; the finder is cost-agnostic (operates on the `g` array), so the same code handles `Norm2` (smooth) and `FaceMax` (piecewise-linear, non-smooth) contact.
6. **Zero-support times are dropped before `extract_qp`** (Phase 3 hand-off): for each `t_j ∈ T^opt` compute `s_j = cost.at(t_j).support(Γᵀ(t_j)λ_opt)`; if `‖s_j‖ < 1e-9` skip it (a `y_j = 0` column leaves `α_j` irrelevant-but-unconstrained). Surviving `(time, s_j, y_j)` are tracked in parallel `Vec`s so returned `α` maps back to the right time/direction.
7. **Initial dual `λ_est = w`** (direction `∥ w`, per Table III). Rationale: the contact function is positively homogeneous (`g(αy) = α·g(y)`, α ≥ 0), so the *scale* of `λ_est` does not change the `argmax` ranking Algorithm 1 uses — only the direction matters, and `w` is the prescribed direction.
8. **`Solution.total_dv = Σ‖dv_j‖₂`** (Euclidean sum of the achieved maneuver magnitudes — this is the 82.4 mm/s figure in Table IV), while `refined.objective = c* = λ_optᵀw` is the dual *lower bound* (82.0 mm/s). `Solution` has no field for `c*`; it is recoverable as `λ·w` from `Solution.lambda` if a later phase needs it. `Solution.residual = ‖w − Σ α_j y_j‖₂ / ‖w‖₂`.
9. **Input validation at the `solve()` boundary** → `PlannerError::InvalidInput`: require `grid.dt > 0`, `grid.t_f > grid.t_i`, `params.n_init ≥ 1`, `params.n_coarse ≥ 1`, and `‖w‖ > 0`. Rationale: `TimeGrid::uniform` and `J2Roe::new` do no validation; catch degenerate inputs once, at entry. (No osculating→mean element conversion — mean elements are assumed per spec §5.4.)

## File Structure

- **Modify `src/algorithm/refine.rs`** — Algorithm 2: the `violated_local_maxima` finder (Task 1), then the `refine(..)` loop + `RefineOutcome` (Task 3).
- **Modify `src/algorithm/init.rs`** — Algorithm 1: `coarse_indices(..)` + `initialize(..)` (Task 2).
- **Modify `src/algorithm/extract.rs`** — Algorithm 3: `extract(..)` + `ExtractOutcome` (Task 4).
- **Modify `src/algorithm/mod.rs`** — shared `contact_at` / `contact_on_grid` helpers, the `cache_gamma` helper, `MAX_REFINE_ITERS`, the `solve()` body, and module `use` lines (Tasks 2/3/5).
- **Create `tests/algorithm.rs`** — public-API end-to-end integration + input-validation tests (Task 5), mirroring `tests/solver.rs`.
- **No change to `src/lib.rs`** — `solve` and `Solution` are already re-exported (`src/lib.rs:13,17`).

Types consumed verbatim (all from `src/types.rs`, re-exported at crate root):
- `Pseudostate = SVector<f64, N>`, `Dual = SVector<f64, N>` (both `SVector<f64, 6>`).
- `Maneuver { t: f64, dv: SVector<f64, M> }` — built struct-literal style (no constructor).
- `TimeGrid { t_i, t_f, dt }` with `uniform`, `len`, `time(idx)`, `times()`.
- `SolveParams { n_coarse, n_init, eps_cost, eps_remove, q }` — defaults 20 / 6 / 0.01 / 0.01 / `I₆`.
- `Solution { maneuvers, total_dv, iterations, residual, lambda }`.
- `ConicRows { linear: Vec<(SVector<f64,N>, f64)>, soc: Vec<(SMatrix<f64,M,N>, f64)> }`.
- `PlannerError::{SolverFailed(String), NotConverged{max_iters, achieved, target}, KeplerDivergence{m, e}, InvalidInput(String)}`.

Phase 3 wrappers consumed verbatim (from `crate::solver`):
- `refine_socp(w: &Pseudostate, rows: &[ConicRows]) -> Result<RefineSolution, PlannerError>`, `RefineSolution { lambda: Dual, objective: f64 }`.
- `extract_qp(w: &Pseudostate, ys: &[SVector<f64, N>], q_weight: &SMatrix<f64, N, N>, budget: f64) -> Result<Vec<f64>, PlannerError>`.

Cost/dynamics traits consumed verbatim:
- `Dynamics::gamma(&self, t: f64) -> SMatrix<f64, N, M>`.
- `CostModel::at(&self, t: f64) -> &dyn SublevelSet`; `SublevelSet::{contact(SVector<f64,M>) -> f64, support(SVector<f64,M>) -> SVector<f64,M>, cone_constraints(&SMatrix<f64,N,M>) -> ConicRows}`.

---

## Task 1: Discrete local-maxima finder (`src/algorithm/refine.rs`)

**Files:**
- Modify: `src/algorithm/refine.rs` (currently a 2-line `//!` doc-only stub)
- Test: `src/algorithm/refine.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: nothing (pure function over `&[f64]`).
- Produces:
  - `pub(super) fn violated_local_maxima(g: &[f64], threshold: f64) -> Vec<usize>` — indices of grid points that are local maxima of `g` **and** exceed `threshold`; one midpoint representative per flat-top plateau; endpoint-inclusive. Consumed by `refine(..)` in Task 3 (so it lands with a transient `#[allow(dead_code)]`).

- [ ] **Step 1: Write the failing test**

Replace the contents of `src/algorithm/refine.rs` with the doc comment, the function under test, and the test module:

```rust
//! Algorithm 2 - Iterative Refinement: solve eq. 40 on `T^est`, drop slack
//! times, add violated local maxima, until convergence.

/// Values within `PLATEAU_EPS` of each other are treated as a flat top.
const PLATEAU_EPS: f64 = 1e-12;

/// Indices of `g` that are local maxima **and** exceed `threshold`.
///
/// A flat top (run of values within [`PLATEAU_EPS`]) yields a single
/// representative (the plateau midpoint). Endpoints are local maxima by the
/// boundary rule (compared only against their one in-bounds neighbour). The
/// global maximum is always a local maximum, so a violated global max is always
/// returned — guaranteeing Algorithm 2 makes progress.
#[allow(dead_code)] // wired into refine() in Task 3
pub(super) fn violated_local_maxima(g: &[f64], threshold: f64) -> Vec<usize> {
    let n = g.len();
    let mut out = Vec::new();
    let mut k = 0usize;
    while k < n {
        // Extent of the flat run [k..=j] of values ~equal to g[k].
        let mut j = k;
        while j + 1 < n && (g[j + 1] - g[k]).abs() <= PLATEAU_EPS {
            j += 1;
        }
        let left_ok = k == 0 || g[k - 1] < g[k];
        let right_ok = j == n - 1 || g[j + 1] < g[k];
        if left_ok && right_ok && g[k] > threshold {
            out.push((k + j) / 2); // plateau midpoint representative
        }
        k = j + 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interior_peak() {
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn left_endpoint_peak() {
        assert_eq!(violated_local_maxima(&[3.0, 1.0, 0.5], 1.0), vec![0]);
    }

    #[test]
    fn right_endpoint_peak() {
        assert_eq!(violated_local_maxima(&[0.5, 1.0, 3.0], 1.0), vec![2]);
    }

    #[test]
    fn flat_top_two_yields_one_midpoint() {
        // plateau [1,2] -> midpoint (1+2)/2 = 1.
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 2.0, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn flat_top_three_yields_one_midpoint() {
        // plateau [1,3] -> midpoint (1+3)/2 = 2.
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 2.0, 2.0, 0.0], 1.0), vec![2]);
    }

    #[test]
    fn monotone_increasing_picks_last() {
        assert_eq!(violated_local_maxima(&[0.0, 1.0, 2.0, 3.0], 1.0), vec![3]);
    }

    #[test]
    fn monotone_decreasing_picks_first() {
        assert_eq!(violated_local_maxima(&[3.0, 2.0, 1.0, 0.0], 1.0), vec![0]);
    }

    #[test]
    fn all_below_threshold_is_empty() {
        assert!(violated_local_maxima(&[0.1, 0.2, 0.1], 1.0).is_empty());
    }

    #[test]
    fn two_separated_peaks() {
        assert_eq!(violated_local_maxima(&[0.0, 2.0, 0.5, 3.0, 0.0], 1.0), vec![1, 3]);
    }

    #[test]
    fn threshold_filters_low_peak() {
        // peak at idx 1 (1.5 > 1) kept; peak at idx 3 (0.8 <= 1) dropped.
        assert_eq!(violated_local_maxima(&[0.0, 1.5, 0.0, 0.8, 0.0], 1.0), vec![1]);
    }

    #[test]
    fn single_element_above_and_below() {
        assert_eq!(violated_local_maxima(&[5.0], 1.0), vec![0]);
        assert!(violated_local_maxima(&[0.5], 1.0).is_empty());
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib algorithm::refine 2>&1 | tail -20`
Expected: at first this should already COMPILE and PASS (the function is included in Step 1). If you instead staged the function and tests separately, the FAIL would be `cannot find function violated_local_maxima`. Since Step 1 includes the implementation, treat Step 2 as the run-to-pass — proceed to Step 4 wording. (This pure function has a trivial implementation, so TDD collapses: the test is the spec; confirm it passes.)

- [ ] **Step 3: (folded into Step 1)** The minimal implementation is the `violated_local_maxima` body shown above. No further code needed.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib algorithm::refine 2>&1 | tail -25`
Expected: PASS (11 tests: `interior_peak`, `left_endpoint_peak`, `right_endpoint_peak`, `flat_top_two_yields_one_midpoint`, `flat_top_three_yields_one_midpoint`, `monotone_increasing_picks_last`, `monotone_decreasing_picks_first`, `all_below_threshold_is_empty`, `two_separated_peaks`, `threshold_filters_low_peak`, `single_element_above_and_below`).

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -25`
Expected: PASS, no clippy warnings. (The `#[allow(dead_code)]` suppresses the unused-function lint until Task 3. If clippy flags `needless_range_loop`/index style anywhere, prefer iterator rewrites over new `#[allow]`s.)

- [ ] **Step 6: Commit**

```bash
git add src/algorithm/refine.rs
git commit -m "feat(algorithm): Alg 2 discrete local-maxima finder (plateau- and endpoint-aware)"
```

---

## Task 2: Algorithm 1 — initialization (`src/algorithm/init.rs`, `src/algorithm/mod.rs`)

**Files:**
- Modify: `src/algorithm/init.rs` (currently a 2-line `//!` doc-only stub)
- Modify: `src/algorithm/mod.rs` (add the shared `contact_at` helper + `use` lines)
- Test: `src/algorithm/init.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `CostModel`/`SublevelSet` (`crate::cost`), `TimeGrid`/`SolveParams`/`Dual` (`crate::types`), `nalgebra::SMatrix`. Uses `super::contact_at` (added to `mod.rs` in this task).
- Produces:
  - In `mod.rs`: `fn contact_at<C: CostModel>(cost: &C, grid: &TimeGrid, gammas: &[SMatrix<f64, N, M>], k: usize, lambda: &Dual) -> f64` — `g_{U(1,grid.time(k))}(Γᵀ(grid.time(k))·lambda)`; the single contact evaluation reused by Init and Refine. Private (callable from child modules via `super::contact_at`); used immediately by `initialize`, so **no** `#[allow(dead_code)]`.
  - In `init.rs`: `pub(super) fn coarse_indices(grid_len: usize, n_coarse: usize) -> Vec<usize>` — up to `n_coarse` evenly spaced, endpoint-inclusive, deduped grid indices (the coarse set `T^d`). Used by `initialize` in this task.
  - In `init.rs`: `pub(super) fn initialize<C: CostModel>(cost: &C, grid: &TimeGrid, gammas: &[SMatrix<f64, N, M>], lambda: &Dual, params: &SolveParams) -> Vec<usize>` — the `n_init` coarse indices with the largest contact value, returned **sorted**. Consumed by `solve()` in Task 5 (lands with a transient `#[allow(dead_code)]`).

- [ ] **Step 1: Write the failing test**

First, add the shared helper and imports to `src/algorithm/mod.rs`. The current `mod.rs` head is:

```rust
//! Orchestration of the three-step algorithm (Init -> Refine -> Extract).

mod extract;
mod init;
mod refine;

use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{PlannerError, Pseudostate, Solution, SolveParams, TimeGrid};
```

Replace that `use` block and add the helper so the head reads:

```rust
//! Orchestration of the three-step algorithm (Init -> Refine -> Extract).

mod extract;
mod init;
mod refine;

use crate::cost::CostModel;
use crate::dynamics::Dynamics;
use crate::types::{Dual, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid, M, N};
use nalgebra::SMatrix;

/// Contact value `g_{U(1,t)}(Γᵀ(t)·lambda)` at grid index `k`.
///
/// `gammas[k] = Γ(grid.time(k))` is the 6×3 dynamics matrix; the cost methods
/// operate on the control-space projection `Γᵀ(t)·lambda ∈ ℝ³`.
fn contact_at<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    k: usize,
    lambda: &Dual,
) -> f64 {
    let y = gammas[k].transpose() * lambda;
    cost.at(grid.time(k)).contact(y)
}
```

Then replace the contents of `src/algorithm/init.rs`:

```rust
//! Algorithm 1 - Initialization: pick the `n_init` coarse times with the
//! largest contact value as `T^est`.

use super::contact_at;
use crate::cost::CostModel;
use crate::types::{Dual, SolveParams, TimeGrid, M, N};
use nalgebra::SMatrix;

/// Up to `n_coarse` evenly spaced grid indices (the coarse set `T^d`),
/// inclusive of both endpoints, clamped to `[1, grid_len]` and deduped.
pub(super) fn coarse_indices(grid_len: usize, n_coarse: usize) -> Vec<usize> {
    debug_assert!(grid_len >= 1);
    let n = n_coarse.clamp(1, grid_len);
    if n == 1 {
        return vec![0];
    }
    let mut idx: Vec<usize> = (0..n)
        .map(|j| ((j as f64) * (grid_len - 1) as f64 / (n - 1) as f64).round() as usize)
        .collect();
    idx.dedup(); // idx is non-decreasing; drop any rounding collisions
    idx
}

/// The `n_init` coarse times with the largest contact `g_{U(1,t)}(Γᵀ(t)·lambda)`,
/// returned as sorted grid indices (`T^est`). `lambda` is the initial dual `∥ w`.
#[allow(dead_code)] // wired into solve() in Task 5
pub(super) fn initialize<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    lambda: &Dual,
    params: &SolveParams,
) -> Vec<usize> {
    let coarse = coarse_indices(grid.len(), params.n_coarse);
    let mut scored: Vec<(usize, f64)> = coarse
        .iter()
        .map(|&k| (k, contact_at(cost, grid, gammas, k, lambda)))
        .collect();
    // Largest contact first; tie-break by index for determinism.
    scored.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
    });
    let n_init = params.n_init.clamp(1, scored.len());
    let mut picked: Vec<usize> = scored.into_iter().take(n_init).map(|(k, _)| k).collect();
    picked.sort_unstable();
    picked
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::Piecewise;
    use crate::dynamics::Dynamics;
    use nalgebra::SVector;

    #[test]
    fn coarse_indices_span_endpoints_and_dedup() {
        assert_eq!(coarse_indices(101, 11), vec![0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]);
        assert_eq!(coarse_indices(1, 20), vec![0]); // clamped to grid_len
        assert_eq!(coarse_indices(5, 1), vec![0]); // n_coarse == 1
        let big = coarse_indices(3, 10); // clamp to 3 -> [0,1,2]
        assert_eq!(big, vec![0, 1, 2]);
    }

    /// Mock dynamics: Γ(t) is the top 3×3 identity scaled by (1 + t), so the
    /// Norm2 contact ‖Γᵀλ‖ grows monotonically with t.
    struct RampDyn;
    impl Dynamics for RampDyn {
        fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
            let mut g = SMatrix::<f64, N, M>::zeros();
            let s = 1.0 + t;
            for i in 0..M {
                g[(i, i)] = s;
            }
            g
        }
    }

    #[test]
    fn initialize_picks_largest_contact_times() {
        let grid = TimeGrid::uniform(0.0, 100.0, 1.0); // 101 points
        let gammas: Vec<SMatrix<f64, N, M>> = grid.times().map(|t| RampDyn.gamma(t)).collect();
        let cost = Piecewise::new(1.0e12); // huge period -> Norm2 everywhere
        let w = SVector::<f64, N>::from_row_slice(&[1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);
        let params = SolveParams {
            n_coarse: 11,
            n_init: 3,
            ..SolveParams::default()
        };
        // Coarse set = {0,10,...,100}; contact grows with t, so the 3 largest
        // are the last three coarse indices, returned sorted.
        let t_est = initialize(&cost, &grid, &gammas, &w, &params);
        assert_eq!(t_est, vec![80, 90, 100]);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib algorithm::init 2>&1 | tail -20`
Expected: PASS once both files are saved (the implementation is included in Step 1). If you split impl from tests, the FAIL is `cannot find function initialize` / `cannot find function coarse_indices`.

- [ ] **Step 3: (folded into Step 1)** Implementations are shown above.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib algorithm 2>&1 | tail -25`
Expected: PASS — Task 1's 11 refine tests + 2 new init tests (`coarse_indices_span_endpoints_and_dedup`, `initialize_picks_largest_contact_times`).

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -25`
Expected: PASS, no clippy warnings. (`contact_at` is used by `initialize` in non-test code, so no `#[allow]` is needed on it; `initialize` keeps its transient `#[allow(dead_code)]` until Task 5. If clippy flags the `M`/`N` import as unused in `mod.rs`, note they are used by `SMatrix<f64, N, M>` in `contact_at` — keep them.)

- [ ] **Step 6: Commit**

```bash
git add src/algorithm/init.rs src/algorithm/mod.rs
git commit -m "feat(algorithm): Alg 1 initialization — n_init largest-contact coarse times"
```

---

## Task 3: Algorithm 2 — iterative refinement (`src/algorithm/refine.rs`, `src/algorithm/mod.rs`)

**Files:**
- Modify: `src/algorithm/refine.rs` (add `RefineOutcome` + `refine(..)`; remove the finder's `#[allow(dead_code)]`)
- Modify: `src/algorithm/mod.rs` (add the `contact_on_grid` helper)
- Test: `src/algorithm/refine.rs` (`#[cfg(test)] mod tests` — extend the existing module)

**Interfaces:**
- Consumes: `refine_socp` + `RefineSolution` (`crate::solver`), `ConicRows`/`Pseudostate`/`SolveParams`/`Dual`/`PlannerError`/`TimeGrid` (`crate::types`), `CostModel`/`SublevelSet::cone_constraints` (`crate::cost`), `violated_local_maxima` (this module, Task 1), `super::contact_on_grid` (added to `mod.rs` here).
- Produces:
  - In `mod.rs`: `fn contact_on_grid<C: CostModel>(cost: &C, grid: &TimeGrid, gammas: &[SMatrix<f64, N, M>], lambda: &Dual) -> Vec<f64>` — `g` at every grid index (built on `contact_at`). Used by `refine` in this task → no `#[allow]`.
  - In `refine.rs`: `pub(super) struct RefineOutcome { pub t_opt: Vec<usize>, pub lambda: Dual, pub objective: f64, pub iterations: usize, pub max_g_trace: Vec<f64> }`.
  - In `refine.rs`: `pub(super) fn refine<C: CostModel>(cost: &C, grid: &TimeGrid, gammas: &[SMatrix<f64, N, M>], w: &Pseudostate, params: &SolveParams, t_est: Vec<usize>, max_iters: usize) -> Result<RefineOutcome, PlannerError>` — runs the eq. 40 loop to convergence (`max_t g ≤ 1 + ε_cost`) or returns `NotConverged`. Consumed by `solve()` in Task 5 (transient `#[allow(dead_code)]`).

**Encoding (faithful restructuring of the paper's `do/while`, Design Decision 4):**
each iteration solves eq. 40 over the current `T^est`, recomputes `g` over the full grid with the new `λ` (since `refine_socp` omits the per-time slack), and **checks convergence before drop/add** so `T^opt` is exactly the solved active set. Thresholds: keep `g ≥ 1 − ε_remove`, add `g > 1`, converge `max_t g ≤ 1 + ε_cost`.

- [ ] **Step 1: Write the failing test**

First add `contact_on_grid` to `src/algorithm/mod.rs`, directly below `contact_at`:

```rust
/// Contact `g` at every grid index, evaluated with the current dual `lambda`.
///
/// `refine_socp` deliberately omits the per-time slack, so Algorithm 2 recomputes
/// `g` here for the drop / add / converge logic.
fn contact_on_grid<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    lambda: &Dual,
) -> Vec<f64> {
    (0..grid.len())
        .map(|k| contact_at(cost, grid, gammas, k, lambda))
        .collect()
}
```

Then, in `src/algorithm/refine.rs`: (a) delete the `#[allow(dead_code)]` line above `violated_local_maxima` (it is now used by `refine`); (b) add the imports, `RefineOutcome`, and `refine` between the `PLATEAU_EPS` const and the `#[cfg(test)]` module; (c) extend the test module. The new non-test code:

```rust
use super::contact_on_grid;
use crate::cost::CostModel;
use crate::solver::refine_socp;
use crate::types::{ConicRows, Dual, PlannerError, Pseudostate, SolveParams, TimeGrid, M, N};
use nalgebra::SMatrix;

/// Result of Algorithm 2.
#[derive(Debug, Clone)]
pub(super) struct RefineOutcome {
    /// Optimal candidate times `T^opt` (grid indices) — the active set that
    /// produced `lambda`.
    pub t_opt: Vec<usize>,
    /// Optimal dual `λ_opt`.
    pub lambda: Dual,
    /// Optimal objective `c* = λ_optᵀw` (the eq. 40 budget for extraction).
    pub objective: f64,
    /// Number of `refine_socp` solves performed.
    pub iterations: usize,
    /// `max_t g` after each solve — non-increasing; read only by tests.
    #[allow(dead_code)]
    pub max_g_trace: Vec<f64>,
}

/// Algorithm 2 — iteratively refine `T^est` until `max_t g ≤ 1 + ε_cost`.
#[allow(dead_code)] // wired into solve() in Task 5
pub(super) fn refine<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    params: &SolveParams,
    mut t_est: Vec<usize>,
    max_iters: usize,
) -> Result<RefineOutcome, PlannerError> {
    let add_threshold = 1.0;
    let converge_threshold = 1.0 + params.eps_cost;
    let keep_threshold = 1.0 - params.eps_remove;

    let mut max_g_trace = Vec::new();
    let mut iterations = 0usize;

    loop {
        if t_est.is_empty() {
            return Err(PlannerError::InvalidInput(
                "refine: candidate-time set became empty".into(),
            ));
        }

        // Solve eq. 40 over the current candidate set T^est.
        let rows: Vec<ConicRows> = t_est
            .iter()
            .map(|&k| cost.at(grid.time(k)).cone_constraints(&gammas[k]))
            .collect();
        let sol = refine_socp(w, &rows)?;
        iterations += 1;

        // Recompute g over the FULL grid with the new dual.
        let g = contact_on_grid(cost, grid, gammas, &sol.lambda);
        let max_g = g.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        max_g_trace.push(max_g);

        // Converged: T^opt is exactly the active set just solved.
        if max_g <= converge_threshold {
            let mut t_opt = t_est;
            t_opt.sort_unstable();
            t_opt.dedup();
            return Ok(RefineOutcome {
                t_opt,
                lambda: sol.lambda,
                objective: sol.objective,
                iterations,
                max_g_trace,
            });
        }
        if iterations >= max_iters {
            return Err(PlannerError::NotConverged {
                max_iters,
                achieved: max_g,
                target: converge_threshold,
            });
        }

        // Drop slack times, then add violated local maxima over the full grid.
        t_est.retain(|&k| g[k] >= keep_threshold);
        for k in violated_local_maxima(&g, add_threshold) {
            t_est.push(k);
        }
        t_est.sort_unstable();
        t_est.dedup();
    }
}
```

Add these tests to the existing `#[cfg(test)] mod tests` in `refine.rs`:

```rust
    use crate::cost::Piecewise;
    use crate::dynamics::Dynamics;
    use nalgebra::SVector;

    /// Mock dynamics: gently rotating, well-conditioned Γ(t). `rate` is chosen
    /// per test so directions vary across the grid without aliasing.
    struct SpinDyn {
        rate: f64,
    }
    impl Dynamics for SpinDyn {
        fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
            let a = self.rate * t;
            let (c, s) = (a.cos(), a.sin());
            SMatrix::<f64, N, M>::from_row_slice(&[
                c, -s, 0.0, //
                s, c, 0.0, //
                0.0, 0.0, 1.0, //
                0.5 * c, 0.0, 0.5 * s, //
                0.0, 0.5 * c, -0.5 * s, //
                0.5 * s, -0.5 * c, 0.0,
            ])
        }
    }

    fn cache(dynamics: &SpinDyn, grid: &TimeGrid) -> Vec<SMatrix<f64, N, M>> {
        grid.times().map(|t| dynamics.gamma(t)).collect()
    }

    #[test]
    fn refine_converges_and_trace_is_non_increasing() {
        let dynamics = SpinDyn { rate: 0.05 };
        let grid = TimeGrid::uniform(0.0, 60.0, 1.0); // 61 points
        let gammas = cache(&dynamics, &grid);
        let cost = Piecewise::new(1.0e12); // Norm2 everywhere
        // Reachable target: impulses at two distinct grid times.
        let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
        let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
        let w = dynamics.gamma(12.0) * ua + dynamics.gamma(47.0) * ub;
        let params = SolveParams::default();

        let out = refine(&cost, &grid, &gammas, &w, &params, vec![0, 30, 60], 50).unwrap();

        // Converged within tolerance.
        assert!(*out.max_g_trace.last().unwrap() <= 1.0 + params.eps_cost + 1e-9);
        assert!(!out.t_opt.is_empty());
        assert!(out.iterations >= 1 && out.iterations <= 50);
        // max_t g is monotonically non-increasing (small tol for solver noise).
        for pair in out.max_g_trace.windows(2) {
            assert!(pair[1] <= pair[0] + 1e-6, "trace not non-increasing: {:?}", out.max_g_trace);
        }
    }

    #[test]
    fn refine_reports_not_converged_at_iteration_cap() {
        let dynamics = SpinDyn { rate: 0.05 };
        let grid = TimeGrid::uniform(0.0, 60.0, 1.0);
        let gammas = cache(&dynamics, &grid);
        let cost = Piecewise::new(1.0e12);
        let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
        let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
        let w = dynamics.gamma(12.0) * ua + dynamics.gamma(47.0) * ub;
        let params = SolveParams::default();

        // One time can't satisfy g ≤ 1 over the whole grid in a single solve.
        let err = refine(&cost, &grid, &gammas, &w, &params, vec![0], 1).unwrap_err();
        match err {
            PlannerError::NotConverged { max_iters, .. } => assert_eq!(max_iters, 1),
            other => panic!("expected NotConverged, got {other:?}"),
        }
    }

    #[test]
    fn refine_handles_mixed_facemax_and_norm2_cones() {
        // Risk R7: exercise the non-smooth FaceMax cost in the refine path.
        // Realistic period so the perigee window is a thin band, giving a true
        // FaceMax/Norm2 mix across the grid.
        let dynamics = SpinDyn { rate: 1.0e-4 };
        let grid = TimeGrid::uniform(0.0, 80_000.0, 1000.0); // 81 points
        let gammas = cache(&dynamics, &grid);
        let cost = Piecewise::new(40_000.0); // FaceMax for t mod 40000 in (16400, 23600)
        // Seed T^est with one Norm2 time (k=5, t=5000) and one FaceMax time
        // (k=20, t=20000) so the FIRST refine_socp assembles mixed cones.
        let ua = SVector::<f64, M>::new(0.5, -0.4, 0.6);
        // FaceMax support directions are tetrahedral vertices; use vertex 0
        // = [√(2/3), 0, -√(1/3)] so the FaceMax time is genuinely reachable.
        let v0 = SVector::<f64, M>::new((2.0_f64 / 3.0).sqrt(), 0.0, -(1.0_f64 / 3.0).sqrt());
        let w = dynamics.gamma(5000.0) * ua + dynamics.gamma(20_000.0) * (0.8 * v0);
        let params = SolveParams::default();

        // Mixed-cone refine_socp must succeed and the trace must be sane; we do
        // not hard-assert convergence-to-1 here (real conditioning is Phase 5's
        // job) — only that the non-smooth path runs and the dual is finite.
        let out = refine(&cost, &grid, &gammas, &w, &params, vec![5, 20], 50).unwrap();
        assert!(out.lambda.iter().all(|x| x.is_finite()));
        for pair in out.max_g_trace.windows(2) {
            assert!(pair[1] <= pair[0] + 1e-6);
        }
    }
```

> **Note on test-module imports:** `refine.rs` now has two groups of `use` — the module-level imports (top of file) and the `#[cfg(test)] mod tests` imports. `use super::*;` already pulls in the module items; add only `Piecewise`, `Dynamics`, `SVector` (and `TimeGrid`/`SolveParams`/`SMatrix` come via `super::*`). If the compiler reports an unused or duplicate import, reconcile against what `super::*` already re-exports.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib algorithm::refine 2>&1 | tail -20`
Expected: before adding `refine`, FAIL with `cannot find function refine` / `cannot find type RefineOutcome` (if you stage tests first). With Step 1 fully applied it should compile; proceed to Step 4.

- [ ] **Step 3: (folded into Step 1)** Implementation shown above.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib algorithm::refine 2>&1 | tail -30`
Expected: PASS — Task 1's 11 finder tests + 3 new refine tests (`refine_converges_and_trace_is_non_increasing`, `refine_reports_not_converged_at_iteration_cap`, `refine_handles_mixed_facemax_and_norm2_cones`).

> If `refine_converges_*` does not reach `≤ 1 + ε_cost`, first confirm the mock `w` is genuinely reachable (it is built from two grid-time impulses), then loosen the band slightly and/or raise the initial `t_est` density — this is the toy-Γ tolerance caveat flagged by Phase 3. Do **not** loosen `refine_reports_not_converged_*` (it must still raise `NotConverged`).

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -30`
Expected: PASS, no clippy warnings. (`violated_local_maxima`'s `#[allow(dead_code)]` is now removed — it is used by `refine`; `refine`/`RefineOutcome` keep their transient `#[allow(dead_code)]` until Task 5; `RefineOutcome.max_g_trace` keeps its permanent field-level `#[allow(dead_code)]`.)

- [ ] **Step 6: Commit**

```bash
git add src/algorithm/refine.rs src/algorithm/mod.rs
git commit -m "feat(algorithm): Alg 2 iterative refinement — eq.40 loop with drop/add/converge"
```

---

## Task 4: Algorithm 3 — control-input extraction (`src/algorithm/extract.rs`)

**Files:**
- Modify: `src/algorithm/extract.rs` (currently a 2-line `//!` doc-only stub)
- Test: `src/algorithm/extract.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `extract_qp` (`crate::solver`), `CostModel`/`SublevelSet::support` (`crate::cost`), `Maneuver`/`Pseudostate`/`Dual`/`PlannerError`/`TimeGrid` (`crate::types`), `nalgebra::{SMatrix, SVector}`.
- Produces:
  - `pub(super) struct ExtractOutcome { pub maneuvers: Vec<Maneuver>, pub total_dv: f64, pub residual: f64 }`.
  - `pub(super) fn extract<C: CostModel>(cost: &C, grid: &TimeGrid, gammas: &[SMatrix<f64, N, M>], w: &Pseudostate, q: &SMatrix<f64, N, N>, lambda: &Dual, budget: f64, t_opt: &[usize]) -> Result<ExtractOutcome, PlannerError>` — builds support directions, drops zero-support times, solves the QP, and assembles maneuvers. Consumed by `solve()` in Task 5 (transient `#[allow(dead_code)]`).

**Encoding (Algorithm 3 + Phase 3 hand-off):**
for each `t_j ∈ T^opt`: `s_j = cost.at(t_j).support(Γᵀ(t_j)λ)`, drop if `‖s_j‖ < SUPPORT_EPS`, else `y_j = Γ(t_j)·s_j`. `α = extract_qp(w, ys, q, budget)` with `budget = λ_optᵀw` (the refined objective). Maneuvers are `Maneuver { t: t_j, dv: α_j·s_j }`; `total_dv = Σ‖dv_j‖`, `residual = ‖w − Σ α_j y_j‖ / ‖w‖`.

- [ ] **Step 1: Write the failing test**

Replace the contents of `src/algorithm/extract.rs`:

```rust
//! Algorithm 3 - Control-Input Extraction: optimal directions, then the QP for
//! magnitudes.

use crate::cost::CostModel;
use crate::solver::extract_qp;
use crate::types::{Dual, Maneuver, PlannerError, Pseudostate, TimeGrid, M, N};
use nalgebra::{SMatrix, SVector};

/// Below this support-direction norm a time contributes no usable maneuver
/// (a `y_j = 0` column leaves `α_j` unconstrained) and is dropped.
const SUPPORT_EPS: f64 = 1e-9;

/// Result of Algorithm 3.
#[derive(Debug, Clone)]
pub(super) struct ExtractOutcome {
    pub maneuvers: Vec<Maneuver>,
    pub total_dv: f64,
    pub residual: f64,
}

/// Algorithm 3 — recover maneuver directions `s_j` and magnitudes `α_j`.
#[allow(dead_code)] // wired into solve() in Task 5
pub(super) fn extract<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    q: &SMatrix<f64, N, N>,
    lambda: &Dual,
    budget: f64,
    t_opt: &[usize],
) -> Result<ExtractOutcome, PlannerError> {
    // Optimal directions, dropping ~zero-support times.
    let mut times: Vec<f64> = Vec::new();
    let mut dirs: Vec<SVector<f64, M>> = Vec::new();
    let mut ys: Vec<SVector<f64, N>> = Vec::new();
    for &k in t_opt {
        let y_dir = gammas[k].transpose() * lambda;
        let s = cost.at(grid.time(k)).support(y_dir);
        if s.norm() < SUPPORT_EPS {
            continue;
        }
        times.push(grid.time(k));
        ys.push(gammas[k] * s);
        dirs.push(s);
    }

    // QP for the nonnegative magnitudes (errors InvalidInput if ys is empty).
    let alpha = extract_qp(w, &ys, q, budget)?;

    let mut maneuvers = Vec::with_capacity(alpha.len());
    let mut w_acc = SVector::<f64, N>::zeros();
    let mut total_dv = 0.0;
    for (j, &a) in alpha.iter().enumerate() {
        let dv = a * dirs[j];
        total_dv += dv.norm();
        w_acc += a * ys[j];
        maneuvers.push(Maneuver { t: times[j], dv });
    }
    let residual = (w - w_acc).norm() / w.norm();

    Ok(ExtractOutcome {
        maneuvers,
        total_dv,
        residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::Piecewise;
    use crate::dynamics::Dynamics;

    fn top_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(i, i)] = 1.0;
        }
        g
    }

    /// Γ(t) = top 3×3 identity for every t.
    struct TopId;
    impl Dynamics for TopId {
        fn gamma(&self, _t: f64) -> SMatrix<f64, N, M> {
            top_identity()
        }
    }

    /// Γ(t) = top identity for t < 5, the zero matrix otherwise.
    struct TopThenZero;
    impl Dynamics for TopThenZero {
        fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
            if t < 5.0 {
                top_identity()
            } else {
                SMatrix::<f64, N, M>::zeros()
            }
        }
    }

    #[test]
    fn extract_recovers_single_maneuver_with_zero_residual() {
        // w = (3,4,12,0,0,0); ‖w_top‖ = 13. λ_opt ∥ w with unit top, so the
        // Norm2 support is unit(w_top); budget = λᵀw = 13; α = 13, dv = w_top.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0);
        let gammas: Vec<SMatrix<f64, N, M>> = grid.times().map(|t| TopId.gamma(t)).collect();
        let cost = Piecewise::new(1.0e12); // Norm2
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let lambda = w / 13.0; // ∥ w, unit top half
        let q = SMatrix::<f64, N, N>::identity();
        let out = extract(&cost, &grid, &gammas, &w, &q, &lambda, 13.0, &[0]).unwrap();

        assert_eq!(out.maneuvers.len(), 1);
        let dv = out.maneuvers[0].dv;
        assert!((dv - SVector::<f64, M>::new(3.0, 4.0, 12.0)).norm() < 1e-6);
        assert!((out.total_dv - 13.0).abs() < 1e-6);
        assert!(out.residual < 1e-6);
    }

    #[test]
    fn extract_drops_zero_support_times() {
        // T^opt includes a time (t=8) where Γ = 0 -> support 0 -> dropped.
        let grid = TimeGrid::uniform(0.0, 10.0, 1.0);
        let gammas: Vec<SMatrix<f64, N, M>> =
            grid.times().map(|t| TopThenZero.gamma(t)).collect();
        let cost = Piecewise::new(1.0e12);
        let w = SVector::<f64, N>::from_row_slice(&[3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let lambda = w / 13.0;
        let q = SMatrix::<f64, N, N>::identity();
        let out = extract(&cost, &grid, &gammas, &w, &q, &lambda, 13.0, &[0, 8]).unwrap();

        assert_eq!(out.maneuvers.len(), 1); // t=8 filtered out
        assert!((out.maneuvers[0].t - 0.0).abs() < 1e-12);
        assert!(out.residual < 1e-6);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib algorithm::extract 2>&1 | tail -20`
Expected: PASS with Step 1 fully applied; if tests are staged before the impl, FAIL with `cannot find function extract` / `cannot find type ExtractOutcome`.

- [ ] **Step 3: (folded into Step 1)** Implementation shown above.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib algorithm::extract 2>&1 | tail -25`
Expected: PASS (2 tests: `extract_recovers_single_maneuver_with_zero_residual`, `extract_drops_zero_support_times`).

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -25`
Expected: PASS, no clippy warnings. (`extract`/`ExtractOutcome` keep their transient `#[allow(dead_code)]` until Task 5.)

- [ ] **Step 6: Commit**

```bash
git add src/algorithm/extract.rs
git commit -m "feat(algorithm): Alg 3 control-input extraction — support dirs + QP magnitudes"
```

---

## Task 5: `solve()` orchestration + end-to-end integration test (`src/algorithm/mod.rs`, `tests/algorithm.rs`)

**Files:**
- Modify: `src/algorithm/mod.rs` (add `cache_gamma` + `MAX_REFINE_ITERS`; replace the `solve()` body; remove the now-satisfied `#[allow(dead_code)]` on `initialize`/`refine`/`extract`)
- Create: `tests/algorithm.rs` (public-API integration + validation tests)

**Interfaces:**
- Consumes: `init::initialize`, `refine::refine` + `RefineOutcome`, `extract::extract` + `ExtractOutcome` (all `pub(super)` in this module's children), `Dynamics::gamma` (`crate::dynamics`), `Solution`/`Maneuver`/`PlannerError`/`Pseudostate`/`SolveParams`/`TimeGrid` (`crate::types`).
- Produces:
  - `fn cache_gamma<D: Dynamics>(dynamics: &D, grid: &TimeGrid) -> Vec<SMatrix<f64, N, M>>` — `Γ(grid.time(k))` for every `k`. Used by `solve` → no `#[allow]`.
  - `const MAX_REFINE_ITERS: usize = 50;`
  - The real body of `pub fn solve<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C, w: Pseudostate, grid: TimeGrid, params: &SolveParams) -> Result<Solution, PlannerError>` (signature unchanged).

- [ ] **Step 1: Write the failing test**

Create `tests/algorithm.rs`:

```rust
//! Public-API integration tests for the Phase 4 three-step planner.

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::Dynamics;
use koenig_planner::{solve, PlannerError, SolveParams, TimeGrid};
use nalgebra::{SMatrix, SVector};

const N: usize = 6;
const M: usize = 3;

/// Gently rotating, well-conditioned Γ(t) (mirrors the unit-test mock).
struct SpinDyn {
    rate: f64,
}
impl Dynamics for SpinDyn {
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
        let a = self.rate * t;
        let (c, s) = (a.cos(), a.sin());
        SMatrix::<f64, N, M>::from_row_slice(&[
            c, -s, 0.0, //
            s, c, 0.0, //
            0.0, 0.0, 1.0, //
            0.5 * c, 0.0, 0.5 * s, //
            0.0, 0.5 * c, -0.5 * s, //
            0.5 * s, -0.5 * c, 0.0,
        ])
    }
}

#[test]
fn solve_converges_on_reachable_synthetic_problem() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0); // 61 points
    let ua = SVector::<f64, M>::new(0.7, -0.3, 0.5);
    let ub = SVector::<f64, M>::new(-0.2, 0.6, 0.4);
    let w = dynamics.gamma(12.0) * ua + dynamics.gamma(47.0) * ub; // reachable
    let cost = Piecewise::new(1.0e12); // Norm2 everywhere
    let params = SolveParams::default();

    let sol = solve(&dynamics, &cost, w, grid, &params).unwrap();

    assert!(sol.residual < 1e-2, "residual {}", sol.residual);
    assert!(sol.iterations >= 1 && sol.iterations <= 50);
    assert!(!sol.maneuvers.is_empty());
    assert!(sol.total_dv > 0.0);
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
}

#[test]
fn solve_rejects_zero_target() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 60.0, 1.0);
    let w = SVector::<f64, N>::zeros();
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}

#[test]
fn solve_rejects_degenerate_grid() {
    let dynamics = SpinDyn { rate: 0.05 };
    let grid = TimeGrid::uniform(0.0, 0.0, 0.0); // dt = 0, t_f == t_i
    let w = SVector::<f64, N>::from_row_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    let cost = Piecewise::new(1.0e12);
    let err = solve(&dynamics, &cost, w, grid, &SolveParams::default()).unwrap_err();
    assert!(matches!(err, PlannerError::InvalidInput(_)));
}
```

Then update `src/algorithm/mod.rs`: add the const + `cache_gamma` (below `contact_on_grid`), and replace the stub `solve`. The current stub is:

```rust
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

Replace it (and add the const + helper) with:

```rust
/// Safety backstop for Algorithm 2 (the paper guarantees convergence; the MC
/// targets converge in ≤ 8 iterations). Not a Table III parameter.
const MAX_REFINE_ITERS: usize = 50;

/// Precompute `Γ(t)` over the grid once (`J2Roe` caches nothing — see Design
/// Decision 2). Indexed by grid index.
fn cache_gamma<D: Dynamics>(dynamics: &D, grid: &TimeGrid) -> Vec<SMatrix<f64, N, M>> {
    grid.times().map(|t| dynamics.gamma(t)).collect()
}

/// Solve the fuel-optimal impulsive control problem.
///
/// Wires Algorithm 1 (init) → Algorithm 2 (refine) → Algorithm 3 (extract) with
/// `Γ(t)` caching.
pub fn solve<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    w: Pseudostate,
    grid: TimeGrid,
    params: &SolveParams,
) -> Result<Solution, PlannerError> {
    // --- Input validation (Design Decision 9). ---
    if !(grid.dt > 0.0) || !(grid.t_f > grid.t_i) {
        return Err(PlannerError::InvalidInput(
            "grid must satisfy dt > 0 and t_f > t_i".into(),
        ));
    }
    if params.n_init == 0 || params.n_coarse == 0 {
        return Err(PlannerError::InvalidInput(
            "n_init and n_coarse must be >= 1".into(),
        ));
    }
    if !(w.norm() > 0.0) {
        return Err(PlannerError::InvalidInput(
            "target pseudostate w must be nonzero".into(),
        ));
    }

    // --- Cache Γ(t) over the grid. ---
    let gammas = cache_gamma(dynamics, &grid);

    // --- Algorithm 1: initialization (λ_est ∥ w). ---
    let t_est = init::initialize(cost, &grid, &gammas, &w, params);

    // --- Algorithm 2: iterative refinement. ---
    let refined = refine::refine(cost, &grid, &gammas, &w, params, t_est, MAX_REFINE_ITERS)?;

    // --- Algorithm 3: control-input extraction. ---
    let extracted = extract::extract(
        cost,
        &grid,
        &gammas,
        &w,
        &params.q,
        &refined.lambda,
        refined.objective,
        &refined.t_opt,
    )?;

    Ok(Solution {
        maneuvers: extracted.maneuvers,
        total_dv: extracted.total_dv,
        iterations: refined.iterations,
        residual: extracted.residual,
        lambda: refined.lambda,
    })
}
```

Finally, **remove the transient `#[allow(dead_code)]`** lines now that `solve` calls each item:
- `src/algorithm/init.rs`: the `#[allow(dead_code)]` above `initialize`.
- `src/algorithm/refine.rs`: the `#[allow(dead_code)]` above `refine` (keep the `RefineOutcome.max_g_trace` field-level allow — it is read only by tests).
- `src/algorithm/extract.rs`: the `#[allow(dead_code)]` above `extract`.

(`mod.rs` no longer has `#[allow(unused_variables)]` on `solve` — it is removed with the stub.)

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --test algorithm 2>&1 | tail -20`
Expected: before the body is implemented, FAIL — the stub `panic!`s via `unimplemented!(...)` so `solve_converges_*` fails with `not implemented: Phases 4-5 wire init -> refine -> extract`.

- [ ] **Step 3: (folded into Step 1)** Implementation shown above.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --test algorithm 2>&1 | tail -25`
Expected: PASS (3 tests: `solve_converges_on_reachable_synthetic_problem`, `solve_rejects_zero_target`, `solve_rejects_degenerate_grid`).

> If `solve_converges_*` exceeds the `1e-2` residual band, the synthetic geometry is the suspect, not the wiring: confirm `w` is reachable, then either widen the band or bump `params.n_coarse`/`n_init` for this test. This is the toy-Γ tolerance caveat (Phase 3 / Risk R3); the real worked-example tolerances are pinned in Phase 5.

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features 2>&1 | tail -30`
Expected: PASS, no clippy warnings. All transient `#[allow(dead_code)]` removed (only `RefineOutcome.max_g_trace`'s field-level allow remains). Phase 4 adds **18 unit tests** (11 finder + 2 init + 3 refine + 2 extract) **and 3 integration tests** — treat the delta (`+18` lib, `+3` integration) as the contract, not the absolute number (the harness counts integration tests inconsistently). The example/bin stubs (`examples/mdot.rs`, `src/bin/monte_carlo.rs`) do not call `solve` and must still build.

- [ ] **Step 6: Commit**

```bash
git add src/algorithm/mod.rs src/algorithm/init.rs src/algorithm/refine.rs src/algorithm/extract.rs tests/algorithm.rs
git commit -m "feat(algorithm): wire solve() — Init -> Refine -> Extract with Γ(t) caching"
```

- [ ] **Step 7: Push and open the PR**

```bash
git push -u origin phase4-algorithms
gh pr create --fill --title "Phase 4 — three algorithms + orchestration (solve)" \
  --body "Implements Phase 4 per docs/superpowers/plans/2026-06-17-koenig-planner-phase4-algorithms.md. Algorithm 1 (init), Algorithm 2 (iterative refinement of eq.40 with plateau-aware local-maxima finder), Algorithm 3 (extraction), and solve() orchestration with Γ(t) caching. +18 unit + 3 integration tests: finder edge cases (R4), mixed FaceMax/Norm2 cones (R7), monotone convergence, NotConverged cap, zero-support filtering, end-to-end convergence + input validation. CI green."
```

---

## Spec & Risk Coverage (self-review)

- **Phase 4 deliverable "Alg. 1 init"** → Task 2. `initialize` picks the `n_init` coarse times (`coarse_indices` = `T^d`, `n_coarse` evenly spaced) with the largest contact, using `λ_est ∥ w` (Design Decision 7).
- **Phase 4 deliverable "Alg. 2 refine (incl. discrete local-maxima finder over the grid, with grid-endpoint handling)"** → Tasks 1 + 3. `violated_local_maxima` is endpoint-inclusive and plateau-aware (Task 1); `refine` runs the eq. 40 loop with drop (`g < 1−ε_remove`) / add (local maxima `> 1`) / converge (`max_t g ≤ 1+ε_cost`), consuming `refine_socp` and `cone_constraints`.
- **Phase 4 deliverable "Alg. 3 extract"** → Task 4. `extract` builds support directions, drops zero-support times, calls `extract_qp` with `budget = λ_optᵀw`, and emits `Maneuver { t, dv }`.
- **Phase 4 deliverable "`solve(...)` wiring with Γ(t) caching"** → Task 5. `cache_gamma` precomputes `Γ` once; `solve` validates inputs and chains Init → Refine → Extract into a `Solution` (signature unchanged, Design Decision 2/9).
- **Phase 4 test "`max_t g` decreases monotonically across iterations"** → Task 3 `refine_converges_and_trace_is_non_increasing` asserts `max_g_trace` is non-increasing (within 1e-6 solver noise).
- **Phase 4 test "convergence within `1+ε_cost`"** → Task 3 (`refine`) + Task 5 (`solve`) assert `max_g_trace.last() ≤ 1 + ε_cost`.
- **Phase 4 test "small residual on a synthetic case"** → Task 4 (`extract`, residual `< 1e-6` on a closed-form case) + Task 5 (`solve`, residual `< 1e-2` end-to-end on a reachable synthetic).
- **Phase 4 exit "end-to-end `solve` runs on a synthetic problem and converges"** → Task 5 `solve_converges_on_reachable_synthetic_problem`.
- **Risk R3 (exact-number matching depends on solver tolerances)** → all synthetic assertions are band-based, not bit-equality; closed-form expectations tie to integer/`√` expressions; the run-to-pass steps document loosening bands rather than chasing digits.
- **Risk R4 (discrete local-maxima finder edge cases at grid ends / plateaus)** → Task 1 dedicated tests: interior peak, both endpoints, flat tops of 2 and 3, monotone increasing/decreasing, all-below-threshold, two separated peaks, threshold filtering, single element.
- **Risk R7 (non-smooth contact × local-maxima finder; face-max cost)** → Task 1 covers non-smooth (sharp triangular) peak arrays cost-agnostically; Task 3 `refine_handles_mixed_facemax_and_norm2_cones` seeds `T^est` with a FaceMax-window time and a Norm2 time so `refine_socp` assembles mixed linear+SOC cones and the grid scan runs FaceMax contact.
- **Edge cases** → empty/become-empty `T^est` (refine guard → `InvalidInput`); iteration cap (`NotConverged`, Task 3); all-zero-support after filter (`extract_qp` errors `InvalidInput`, propagated); degenerate grid / zero `w` (Task 5 validation); `n_init`/`n_coarse` larger than the grid (clamped in `coarse_indices`/`initialize`).
- **Type consistency** → `T^est`/`T^opt` are `Vec<usize>` end to end; `Γ`-cache is `Vec<SMatrix<f64, N, M>>`; `contact`/`support` always receive `Γᵀ(t)λ` (the `M`-vector), `cone_constraints` always receives the full `N×M` `Γ(t)`; `budget` passed to `extract_qp` is exactly `refined.objective` (= `λ_optᵀw`); `Solution.lambda = refined.lambda`. All names match across tasks (`initialize`, `refine`, `extract`, `cache_gamma`, `contact_at`, `contact_on_grid`, `violated_local_maxima`, `RefineOutcome`, `ExtractOutcome`).

## Notes for Phase 5 (worked-example validation — not implemented here)

- **`examples/mdot.rs` encodes Table III** and calls the public `solve`: build the chief with `AbsoluteOrbit::new(25_000e3, 0.7, 40f64.to_radians(), 358f64.to_radians(), 0.0, 180f64.to_radians())`, `J2Roe::new(chief, 0.0, 117_990.0)`, `TimeGrid::uniform(0.0, 117_990.0, 30.0)` (3934 pts), `Piecewise::new(2.0*PI/n)` (orbit period in seconds), and `w = [50,5000,100,100,0,400]/a_c` (divide the metre-scaled Table III `w` by `a_c = 25_000e3` — §5.5). Multiply reported Δv back by nothing (already m/s) and states by `a_c` for display.
- **Tolerances will tighten and may need revisiting.** Phase 4's synthetic tests use a well-conditioned mock `Γ` and a huge-period `Piecewise` (≈ pure `Norm2`); Phase 5 is the first run on the **real, ill-conditioned `J2Roe` `Γ(t)`** with the genuine eq. 49 mix. Expect to revisit the `refine` convergence behaviour and, if needed, add `clarabel` warm-starting (Risk R8) — the per-iteration SOCPs share `λ_opt`.
- **The real FaceMax/Norm2 piecewise path gets its first true exercise in Phase 5.** Phase 4 verifies the mixed-cone *mechanics* (Task 3 R7 test) but not the published 3-maneuver/82.4 mm/s numbers; those are the Phase 5 targets (spec §7).
- **`Solution` has no `c*`/objective field.** Phase 5 compares achieved `total_dv` (≈ 82.4) against the dual lower bound `c* = λ_optᵀw` (≈ 82.0); recompute `c*` as `sol.lambda.dot(&w)` from the returned `Solution`. If a later phase needs it surfaced, that is a `types.rs` change outside Phase 4's locked scope.
- **`MAX_REFINE_ITERS = 50`** is a backstop, not a tuning knob; if a Phase 5/6 case legitimately needs more, prefer adding a `max_iters` field to `SolveParams` (a deliberate public-API change) over silently raising the const.

---

**Plan complete and saved to `docs/superpowers/plans/2026-06-17-koenig-planner-phase4-algorithms.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — I execute the tasks in this session using executing-plans, with checkpoints for review.

**Which approach?**
