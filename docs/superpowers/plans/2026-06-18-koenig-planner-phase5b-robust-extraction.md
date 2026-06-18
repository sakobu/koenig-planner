# Phase 5b — Robust Algorithm-3 Extraction (Direct Min-Fuel SOCP) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Algorithm 3 (control-input extraction) robust on the degenerate flat contacts of e=0.7 orbits by replacing the fixed-support-direction magnitude QP with a direct min-fuel SOCP that recovers full 3-DOF maneuvers, driving the worked example to <0.1% residual and the Hunter L2 case to <0.01%.

**Architecture:** The refinement dual (`refine_socp`, eq. 40) is already correct; only the *primal recovery* is weak. We add a new solver `min_fuel_socp` that solves `min Σⱼ f_{tⱼ}(Δvⱼ) s.t. Σⱼ Γ(tⱼ)·Δvⱼ = w` over the converged active set `T^opt`, with each maneuver charged by its *true per-time cost* (L2 second-order cone for `Norm2` times; a nonnegative-combination LP over the `V_vertex` directions for `FaceMax` times). By conic strong duality this attains the dual value `c*`, and because `Δvⱼ` is a full vector (not a fixed support direction) the active times now span `w`, so the residual collapses. Sum-of-norms is group-sparse, so few maneuvers come out nonzero. The cost layer gains one new trait method (`fuel_generator`) describing how each cost builds and charges a Δv.

**Tech Stack:** Rust 2021, `nalgebra` 0.35 (static-dim `SMatrix`/`SVector`), `clarabel` 0.11 (interior-point conic solver: zero/nonnegative/second-order cones), `thiserror` 2.0.

## Global Constraints

- **Dynamics are FROZEN. Do NOT touch `src/dynamics/**`.** They are independently finite-difference verified (`tests/fd_stm.rs`, `tests/fd_b_matrix.rs`) at the Koenig chief, the e=0.3 fixture, and the Hunter chief. If any dynamics file changes, those two tests are the gate. The `Φ₂₁` dt² typo is already fixed (`src/dynamics/stm.rs`).
- **The paper's Table III/IV / 82.4 mm/s figures are NOT reproducible** (internally inconsistent with the corrected dynamics — see spec §6 Phase 5 and `tests/worked_example.rs::paper_table_iv_does_not_reconstruct`). Do not chase them. Headline numbers are *our* FD-verified optimum (~80.9 mm/s worked example; ~2.48e-4 m/s Hunter L2), validated by self-consistency against the dual.
- **Versions (verbatim):** `nalgebra = "0.35"`, `clarabel = "0.11"`, `thiserror = "2.0"`, `rust-version = "1.92"`. No new dependencies.
- **Dimensions:** `N = 6` (ROE state), `M = 3` (RTN control). Both are `pub const` in `src/types.rs`.
- **clarabel conventions (from Phase 3):** solves `min ½xᵀPx + qᵀx s.t. Ax + s = b, s ∈ K`. `DefaultSolver::new(&P,&q,&A,&b,&cones,settings)`; `solver.solve()` **requires `use clarabel::solver::IPSolver` in scope**; read `solver.solution.{x, status, obj_val}`. `CscMatrix::from(&dense)` (a `&Vec<Vec<f64>>`) drops exact zeros. Accept `SolverStatus::Solved | AlmostSolved` (use the existing `check_status`). Use the existing `silent_settings()`. Cones: `ZeroConeT(n)` (equalities), `NonnegativeConeT(n)`, `SecondOrderConeT(n)` — all from `clarabel::solver::*`. The order of `cones` MUST match the row order of `A`/`b`.
- **CI gate (all `--all-features`):** `cargo fmt --check` + `cargo clippy --all-features -- -D warnings` + `cargo build --all-features` + `cargo test --all-features`. The Linux runner installs `libfontconfig1-dev` for the `plotters` `validation` feature.
- **Lint gotcha:** helpers used only by `#[cfg(test)]` until a consumer lands trip `dead_code` under `clippy -D warnings` — use `#[allow(dead_code)]`, **NOT `#[expect(dead_code)]`** (which mis-fires in the `cfg(test)` build). For non-finite guards use `!x.is_finite()`, never `x <= 0.0` (admits NaN) or `!(x > 0.0)` (trips `clippy::neg_cmp_op_on_partial_ord`).
- **Process:** this repo **requires feature-branch PRs even for docs** (direct pushes to `main` are blocked). Work on branch `phase5b-robust-extraction`.

## Design Decisions (locked)

- **D1 — Faithful, gauge-aware min-fuel (chosen) vs L2-only (rejected).** We charge each maneuver by its *true* per-time cost: `Norm2` → L2 (SOC), `FaceMax` → polytopic gauge (LP over `V_vertex` columns). This keeps the primal cost equal to the dual `c*` *at any optimum* (the dual was computed with the piecewise cost), so the dual/primal self-consistency check that the whole phase rests on stays valid even if an optimal impulse lands at a `FaceMax` time. The simpler "min Σ‖vⱼ‖₂ for all times" form (spec fix #1, literal) is a valid descope that would also pass the two acceptance cases (the worked example's optimum is all-`Norm2`, Hunter is pure L2), but it silently mis-costs any `FaceMax`-active maneuver — rejected for a *faithful* reimplementation.
- **D2 — Candidate set = `T^opt`** (the converged active set returned by `refine`). On the degenerate flat contact `T^opt` is large and rich, so its `Γ` columns span ℝ⁶ ⇒ the equality `Σ Γ Δv = w` is feasible and the min equals `c*`. (The optimal primal measure is supported within `T^opt`, so no enrichment is needed; if a future case shows otherwise, augment with the distinct contact peaks — not done here.)
- **D3 — `extract` signature changes** (it is `pub(super)`, internal): drop `q` (the QP weight — the min-fuel SOCP does not use it) and drop `lambda` (no longer needed: full-DOF recovery does not use support directions). Keep `budget` (= the dual `c*` from `refine`) and use it only as a self-consistency *check*. `SolveParams.q` is **kept** (public API stability) but is now unused by the extractor; document it as reserved.
- **D4 — `extract_qp` is retained** (public, re-exported, 8 closed-form tests) as the faithful Algorithm-3-as-printed QP primitive, but is **no longer on the `solve()` path**. Do not delete it.
- **D5 — Pruning is characterization-first.** Interior-point solutions are not exactly sparse; inactive maneuvers come back at ~solver-tolerance magnitude. Prune by a *relative* threshold and set its exact value (and the "small maneuver set" bound) by observing the real magnitude distribution before locking the test (Phase-5 methodology: run, observe, then assert).

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `src/types.rs` | core value types | **Modify**: add `FuelGenerator` enum |
| `src/cost/mod.rs` | `SublevelSet`/`CostModel` traits | **Modify**: add `fuel_generator` trait method |
| `src/cost/norm2.rs` | L2 cost | **Modify**: implement `fuel_generator` → `Norm` |
| `src/cost/facemax.rs` | FaceMax cost | **Modify**: implement `fuel_generator` → `Polytope`; expose vertex columns |
| `src/solver/min_fuel_socp.rs` | the new direct min-fuel SOCP | **Create** |
| `src/solver/mod.rs` | solver wrappers + helpers | **Modify**: declare + re-export `min_fuel_socp` |
| `src/lib.rs` | public surface | **Modify**: re-export `min_fuel_socp`, `MinFuelSolution`, `FuelGenerator` |
| `src/algorithm/extract.rs` | Algorithm 3 | **Modify**: call `min_fuel_socp`, prune, new signature |
| `src/algorithm/mod.rs` | `solve` orchestration | **Modify**: update the `extract::extract(...)` call |
| `src/algorithm/refine.rs` | Algorithm 2 | **Modify**: add `active_set_trace` field |
| `examples/mdot.rs` | worked-example report | **Modify**: report the now-low residual / small plan |
| `tests/worked_example.rs` | §7 self-consistency gate | **Modify**: tighten residual + add Hunter L2 test |
| `tests/algorithm.rs` | public-API integration tests | **Modify**: add the real-`J2Roe` drop-then-readd refine test |
| `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md` | spec | **Modify**: mark Phase 5b ✅ Done |

---

## Task 1: `FuelGenerator` type + `SublevelSet::fuel_generator`

Add the descriptor the min-fuel SOCP needs from each cost, and implement it for both sublevel sets. The trait method addition forces both impls to change in the same commit (so the crate keeps compiling).

**Files:**
- Modify: `src/types.rs` (add `FuelGenerator`, after `ConicRows`)
- Modify: `src/cost/mod.rs` (add trait method + import)
- Modify: `src/cost/norm2.rs` (impl)
- Modify: `src/cost/facemax.rs` (impl + expose vertex columns)

**Interfaces:**
- Produces: `enum FuelGenerator { Norm, Polytope(Vec<SVector<f64, M>>) }` (in `crate::types`, re-exported from `crate`); trait method `fn fuel_generator(&self) -> FuelGenerator` on `SublevelSet`.
- Consumes: nothing (leaf task).

- [ ] **Step 1: Write the failing tests**

In `src/cost/norm2.rs`, inside `mod tests`, add:

```rust
    #[test]
    fn fuel_generator_is_norm() {
        use crate::types::FuelGenerator;
        assert_eq!(Norm2.fuel_generator(), FuelGenerator::Norm);
    }
```

In `src/cost/facemax.rs`, inside `mod tests`, add:

```rust
    #[test]
    fn fuel_generator_is_polytope_of_unit_vertices() {
        use crate::types::FuelGenerator;
        match FaceMax.fuel_generator() {
            FuelGenerator::Polytope(dirs) => {
                assert_eq!(dirs.len(), 4);
                // Same four unit tetrahedral directions used by contact/support.
                for (d, v) in dirs.iter().zip(vertex_columns()) {
                    assert_relative_eq!(*d, v, epsilon = 1e-12);
                    assert_relative_eq!(d.norm(), 1.0, epsilon = 1e-12);
                }
            }
            other => panic!("expected Polytope, got {other:?}"),
        }
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --lib cost:: 2>&1 | tail -20`
Expected: FAIL — `no method named fuel_generator` / `no variant FuelGenerator`.

- [ ] **Step 3: Add the `FuelGenerator` type**

In `src/types.rs`, immediately after the `ConicRows` struct (around line 122), add:

```rust
/// Primal fuel generator for one maneuver in the direct min-fuel SOCP
/// (Phase 5b, Algorithm 3). Describes how a Δv at one candidate time is built
/// from solver variables and how it is charged, mirroring the cost's unit
/// sublevel set:
///
/// * `Norm` — a free vector `v ∈ ℝᴹ` charged its L2 norm `‖v‖₂` (an `‖v‖ ≤ c`
///   second-order cone). This is the L2 cost.
/// * `Polytope(dirs)` — a nonnegative combination `Δv = Σₖ θₖ·dirs[k]`,
///   `θₖ ≥ 0`, charged `Σₖ θₖ` (a nonnegative-cone LP). The unit ball is
///   `conv{0, dirs…}`, so this is the gauge of that polytope; `dirs` are the
///   FaceMax `V_vertex` columns.
#[derive(Debug, Clone, PartialEq)]
pub enum FuelGenerator {
    /// L2: free `v ∈ ℝᴹ` charged `‖v‖₂`.
    Norm,
    /// Polytopic: `Δv = Σₖ θₖ·dirs[k]`, `θ ≥ 0`, charged `Σₖ θₖ`.
    Polytope(Vec<SVector<f64, M>>),
}
```

- [ ] **Step 4: Add the trait method**

In `src/cost/mod.rs`, add `FuelGenerator` to the import and a method to `SublevelSet`:

```rust
use crate::types::{ConicRows, FuelGenerator, M, N};
```

```rust
pub trait SublevelSet {
    /// Contact function `g(y) = max_{z in U} y . z`.
    fn contact(&self, y: SVector<f64, M>) -> f64;
    /// Support direction `s(y) = argmax_{z in U} y . z`.
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M>;
    /// Conic rows encoding `g_{U(1,t)}(Gamma^T(t) lambda) <= 1`.
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows;
    /// Primal fuel generator for the direct min-fuel SOCP (Phase 5b): how a Δv
    /// in this sublevel set is built from solver variables and charged.
    fn fuel_generator(&self) -> FuelGenerator;
}
```

- [ ] **Step 5: Implement for `Norm2`**

In `src/cost/norm2.rs`, add to `impl SublevelSet for Norm2` (and extend the import):

```rust
use crate::types::{ConicRows, FuelGenerator, M, N};
```

```rust
    fn fuel_generator(&self) -> FuelGenerator {
        FuelGenerator::Norm
    }
```

- [ ] **Step 6: Implement for `FaceMax`**

In `src/cost/facemax.rs`, extend the import and add the impl method. `vertex_columns()` already exists (private module fn); reuse it:

```rust
use crate::types::{ConicRows, FuelGenerator, M, N};
```

```rust
    fn fuel_generator(&self) -> FuelGenerator {
        // Unit ball U(1) = conv{0, V_vertex columns}; its gauge is
        //   f(v) = min{ Σₖ θₖ : Σₖ θₖ vₖ = v, θ ≥ 0 }.
        // (The eq.48 V_face matrix is the inconsistent printed form and is used
        // nowhere; the algorithm's geometry is the V_vertex columns.)
        FuelGenerator::Polytope(vertex_columns().to_vec())
    }
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `cargo test --lib cost:: 2>&1 | tail -20`
Expected: PASS (both new tests, plus all existing cost tests still green).

- [ ] **Step 8: Commit**

```bash
git add src/types.rs src/cost/mod.rs src/cost/norm2.rs src/cost/facemax.rs
git commit -m "feat(cost): add FuelGenerator + SublevelSet::fuel_generator for Phase 5b min-fuel SOCP"
```

---

## Task 2: `min_fuel_socp` solver

The new direct min-fuel SOCP. Closed-form unit tests include the **degenerate collinear-support case** that the magnitude QP cannot solve but this solver can.

**Files:**
- Create: `src/solver/min_fuel_socp.rs`
- Modify: `src/solver/mod.rs` (declare + re-export)
- Modify: `src/lib.rs` (re-export)
- Test: in-module `#[cfg(test)] mod tests` in `src/solver/min_fuel_socp.rs`

**Interfaces:**
- Consumes: `FuelGenerator` (Task 1); `check_status`, `silent_settings` (existing, `pub(crate)` in `src/solver/mod.rs`).
- Produces:
  - `struct MinFuelSolution { pub dvs: Vec<SVector<f64, M>>, pub objective: f64 }`
  - `fn min_fuel_socp(w: &Pseudostate, gammas: &[SMatrix<f64, N, M>], generators: &[FuelGenerator]) -> Result<MinFuelSolution, PlannerError>` — `dvs` aligned to the input candidate times (one per `gammas`/`generators` entry); `objective` = total fuel `Σⱼ f_{tⱼ}(Δvⱼ)`.

- [ ] **Step 1: Write the failing tests**

Create `src/solver/min_fuel_socp.rs` with ONLY the test module first (so the file compiles enough to fail meaningfully). Paste the full tests block from Step 4 below, then run Step 2. (If you prefer, write the implementation and tests together — but run the tests-fail check first by stubbing `min_fuel_socp` to `unimplemented!()`.)

- [ ] **Step 2: Run to verify failure**

Run: `cargo test --lib min_fuel_socp 2>&1 | tail -20`
Expected: FAIL (unresolved `min_fuel_socp` / `MinFuelSolution`, or `unimplemented!`).

- [ ] **Step 3: Write the implementation**

Replace the file contents of `src/solver/min_fuel_socp.rs` with:

```rust
//! Direct min-fuel SOCP (Phase 5b): recover full 3-DOF maneuvers over a fixed
//! candidate-time set by minimizing the (piecewise) fuel cost subject to exact
//! reachability `Σⱼ Γ(tⱼ)·Δvⱼ = w`.
//!
//! This replaces the fixed-support-direction magnitude QP (`extract_qp`) inside
//! Algorithm 3, which is not robust on the degenerate flat contacts of e=0.7
//! orbits: there the per-time support directions are nearly collinear and cannot
//! span `w`. Freeing each maneuver to a full vector (charged by its true cost)
//! recovers `w` to ~0 residual; by conic strong duality the optimum equals the
//! eq.40 dual value `c*`. Sum-of-norms is group-sparse, so few maneuvers come
//! out nonzero.

use crate::solver::{check_status, silent_settings};
use crate::types::{FuelGenerator, PlannerError, Pseudostate, M, N};
use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT, SupportedConeT, ZeroConeT,
};
use nalgebra::{SMatrix, SVector};

/// Result of the direct min-fuel SOCP.
#[derive(Debug, Clone)]
pub struct MinFuelSolution {
    /// Recovered Δv per input candidate time (aligned to `gammas`/`generators`;
    /// times the optimum does not use come back ≈ 0 via the sum-of-norms penalty).
    pub dvs: Vec<SVector<f64, M>>,
    /// Total fuel cost `Σⱼ f_{tⱼ}(Δvⱼ)` (≈ the dual budget `c*`).
    pub objective: f64,
}

/// Solve `min Σⱼ f_{tⱼ}(Δvⱼ) s.t. Σⱼ Γ(tⱼ)·Δvⱼ = w` over the candidate times
/// described pairwise by (`gammas`, `generators`).
///
/// `f_{tⱼ}` is the cost whose unit ball the generator describes: an L2 norm
/// (`Norm` → one second-order cone) or a polytopic gauge (`Polytope` →
/// nonnegative-cone LP over the vertex directions). Returns one Δv per candidate
/// time plus the optimal fuel.
pub fn min_fuel_socp(
    w: &Pseudostate,
    gammas: &[SMatrix<f64, N, M>],
    generators: &[FuelGenerator],
) -> Result<MinFuelSolution, PlannerError> {
    let k = gammas.len();
    if k == 0 || k != generators.len() {
        return Err(PlannerError::InvalidInput(
            "min_fuel_socp: need >= 1 candidate time and matching generators".into(),
        ));
    }

    // --- Per-maneuver variable layout (contiguous block per maneuver). ---
    //   Norm        -> [c_j, v_j(M)]   (M+1 vars; Δv_j = v_j; cost = c_j ≥ ‖v_j‖)
    //   Polytope(p) -> [θ_j(p)]        (p   vars; Δv_j = Σₖ θ_jk d_k; cost = Σ θ)
    let mut offsets = Vec::with_capacity(k);
    let mut n_var = 0usize;
    for gen in generators {
        offsets.push(n_var);
        n_var += match gen {
            FuelGenerator::Norm => M + 1,
            FuelGenerator::Polytope(dirs) => dirs.len(),
        };
    }

    // --- Objective q (P = 0): minimize the sum of the cost variables. ---
    let mut q = vec![0.0f64; n_var];
    for (j, gen) in generators.iter().enumerate() {
        match gen {
            FuelGenerator::Norm => q[offsets[j]] = 1.0, // the epigraph c_j
            FuelGenerator::Polytope(dirs) => {
                for kk in 0..dirs.len() {
                    q[offsets[j] + kk] = 1.0;
                }
            }
        }
    }

    // --- Equality-constraint column vectors per variable: the 6-vector that
    //     multiplies each variable in `Σ Γ Δv = w`. For Norm, the c_j column is
    //     0 and the v_j columns are Γ[:,m]; for Polytope, column k is Γ·d_k. ---
    let mut eq_cols: Vec<Vec<SVector<f64, N>>> = Vec::with_capacity(k);
    for (j, gen) in generators.iter().enumerate() {
        let g = &gammas[j];
        match gen {
            FuelGenerator::Norm => {
                let mut cols = vec![SVector::<f64, N>::zeros()]; // c_j
                for m in 0..M {
                    cols.push(SVector::<f64, N>::from_iterator(g.column(m).iter().copied()));
                }
                eq_cols.push(cols);
            }
            FuelGenerator::Polytope(dirs) => {
                eq_cols.push(dirs.iter().map(|d| g * d).collect());
            }
        }
    }

    // --- Assemble A, b, cones in lockstep:
    //     (1) N equality rows  -> ZeroCone(N)
    //     (2) per maneuver:  Norm -> SOC(M+1) on (c_j, v_j);  Polytope -> Nonneg(p)
    let mut a_rows: Vec<Vec<f64>> = Vec::new();
    let mut b: Vec<f64> = Vec::new();
    let mut cones: Vec<SupportedConeT<f64>> = Vec::new();

    // (1) Equality block: A·x = w  (s ∈ ZeroCone).
    for r in 0..N {
        let mut row = vec![0.0; n_var];
        for j in 0..k {
            for (c, col) in eq_cols[j].iter().enumerate() {
                row[offsets[j] + c] = col[r];
            }
        }
        a_rows.push(row);
        b.push(w[r]);
    }
    cones.push(ZeroConeT(N));

    // (2) Per-maneuver cost cones.  clarabel: s = b - A·x.
    for (j, gen) in generators.iter().enumerate() {
        match gen {
            FuelGenerator::Norm => {
                // s = (c_j, v_j) ∈ SOC(M+1): A = -selector, b = 0.
                let mut s0 = vec![0.0; n_var];
                s0[offsets[j]] = -1.0; // s_0 = c_j
                a_rows.push(s0);
                b.push(0.0);
                for m in 0..M {
                    let mut row = vec![0.0; n_var];
                    row[offsets[j] + 1 + m] = -1.0; // s_{1+m} = v_{j,m}
                    a_rows.push(row);
                    b.push(0.0);
                }
                cones.push(SecondOrderConeT(M + 1));
            }
            FuelGenerator::Polytope(dirs) => {
                // s_k = θ_jk ≥ 0.
                for kk in 0..dirs.len() {
                    let mut row = vec![0.0; n_var];
                    row[offsets[j] + kk] = -1.0;
                    a_rows.push(row);
                    b.push(0.0);
                }
                cones.push(NonnegativeConeT(dirs.len()));
            }
        }
    }

    let a_csc = CscMatrix::from(&a_rows);
    let p_csc = CscMatrix::<f64>::zeros((n_var, n_var));

    let mut solver = DefaultSolver::new(&p_csc, &q, &a_csc, &b, &cones, silent_settings())
        .map_err(|e| PlannerError::SolverFailed(format!("clarabel setup failed: {e:?}")))?;
    solver.solve();
    check_status(solver.solution.status)?;
    let x = &solver.solution.x;

    // --- Reconstruct Δv per maneuver and total fuel. ---
    let mut dvs = Vec::with_capacity(k);
    let mut objective = 0.0;
    for (j, gen) in generators.iter().enumerate() {
        match gen {
            FuelGenerator::Norm => {
                let v = SVector::<f64, M>::from_iterator((0..M).map(|m| x[offsets[j] + 1 + m]));
                objective += x[offsets[j]]; // c_j ≈ ‖v‖
                dvs.push(v);
            }
            FuelGenerator::Polytope(dirs) => {
                let mut v = SVector::<f64, M>::zeros();
                for (kk, d) in dirs.iter().enumerate() {
                    let theta = x[offsets[j] + kk];
                    v += theta * d;
                    objective += theta;
                }
                dvs.push(v);
            }
        }
    }

    Ok(MinFuelSolution { dvs, objective })
}
```

- [ ] **Step 4: Add the test module**

Append to `src/solver/min_fuel_socp.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Γ (6x3) with top 3x3 = I, bottom 3x3 = 0: Γv = (v; 0), Γᵀλ = (λ₁,λ₂,λ₃).
    fn gamma_top_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(i, i)] = 1.0;
        }
        g
    }

    fn norm_gen() -> FuelGenerator {
        FuelGenerator::Norm
    }

    // The four FaceMax V_vertex columns (must match src/cost/facemax.rs).
    fn vertex_dirs() -> Vec<SVector<f64, M>> {
        let a = (2.0_f64 / 3.0).sqrt();
        let b = (1.0_f64 / 3.0).sqrt();
        vec![
            SVector::<f64, M>::new(a, 0.0, -b),
            SVector::<f64, M>::new(-a, 0.0, -b),
            SVector::<f64, M>::new(0.0, a, b),
            SVector::<f64, M>::new(0.0, -a, b),
        ]
    }

    fn w6(v: [f64; N]) -> Pseudostate {
        Pseudostate::from_row_slice(&v)
    }

    fn residual(w: &Pseudostate, gammas: &[SMatrix<f64, N, M>], dvs: &[SVector<f64, M>]) -> f64 {
        let mut acc = SVector::<f64, N>::zeros();
        for (g, dv) in gammas.iter().zip(dvs) {
            acc += g * dv;
        }
        (w - acc).norm() / w.norm()
    }

    #[test]
    fn empty_or_mismatched_input_is_invalid() {
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(matches!(
            min_fuel_socp(&w, &[], &[]).unwrap_err(),
            PlannerError::InvalidInput(_)
        ));
        assert!(matches!(
            min_fuel_socp(&w, &[gamma_top_identity()], &[]).unwrap_err(),
            PlannerError::InvalidInput(_)
        ));
    }

    #[test]
    fn single_norm_time_recovers_exact_dv() {
        // w = (3,4,12,0,0,0), Γ = [I;0]: Δv = (3,4,12), fuel = ‖Δv‖ = 13.
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let gammas = [gamma_top_identity()];
        let sol = min_fuel_socp(&w, &gammas, &[norm_gen()]).unwrap();
        assert_relative_eq!(sol.dvs[0], SVector::<f64, M>::new(3.0, 4.0, 12.0), epsilon = 1e-5);
        assert_relative_eq!(sol.objective, 13.0, epsilon = 1e-4);
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
    }

    #[test]
    fn degenerate_collinear_support_still_spans_w() {
        // THE Phase-5b case in miniature. Two times whose Γ map control space to
        // DIFFERENT pseudostate planes; the per-time support directions of a
        // shared λ are collinear (the failure mode of extract_qp), yet w needs a
        // contribution from each time. Full-DOF min-fuel recovers w exactly.
        // Time A: Γ_A = [I;0] (controls coords 0..2). Time B: Γ_B maps to coords 3..5.
        let mut gb = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            gb[(M + i, i)] = 1.0; // bottom identity
        }
        let gammas = [gamma_top_identity(), gb];
        let w = w6([1.0, 0.0, 0.0, 0.0, 2.0, 0.0]); // needs BOTH times
        let sol = min_fuel_socp(&w, &gammas, &[norm_gen(), norm_gen()]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5, "residual too large");
        // Fuel = ‖(1,0,0)‖ + ‖(0,2,0)‖ = 1 + 2 = 3.
        assert_relative_eq!(sol.objective, 3.0, epsilon = 1e-4);
    }

    #[test]
    fn unused_time_comes_back_zero() {
        // Second time contributes nothing to a reachable w -> its Δv ≈ 0.
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let gammas = [gamma_top_identity(), gamma_top_identity()];
        let sol = min_fuel_socp(&w, &gammas, &[norm_gen(), norm_gen()]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        // Total fuel is still 13 (split across the two identical times any way).
        assert_relative_eq!(sol.objective, 13.0, epsilon = 1e-3);
    }

    #[test]
    fn facemax_single_vertex_costs_unit() {
        // w = Γ·v0 with Γ = [I;0], v0 the first vertex direction (a unit vector).
        // The gauge of a vertex is 1, so fuel = 1 and Δv ≈ v0.
        let dirs = vertex_dirs();
        let gammas = [gamma_top_identity()];
        let mut w = SVector::<f64, N>::zeros();
        for r in 0..M {
            w[r] = dirs[0][r];
        }
        let sol = min_fuel_socp(&w, &gammas, &[FuelGenerator::Polytope(dirs.clone())]).unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        assert_relative_eq!(sol.objective, 1.0, epsilon = 1e-3);
        assert_relative_eq!(sol.dvs[0], dirs[0], epsilon = 1e-3);
    }

    #[test]
    fn mixed_norm_and_facemax_times() {
        // One Norm time (coords 0..2) + one FaceMax time (coords 3..5). Both used.
        let dirs = vertex_dirs();
        let mut gb = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            gb[(M + i, i)] = 1.0;
        }
        let gammas = [gamma_top_identity(), gb];
        // w: (1,0,0) on the Norm side; gb·v2 (a vertex) on the FaceMax side.
        let mut w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let face_part = gb * dirs[2];
        for r in 0..N {
            w[r] += face_part[r];
        }
        let sol = min_fuel_socp(
            &w,
            &gammas,
            &[norm_gen(), FuelGenerator::Polytope(dirs.clone())],
        )
        .unwrap();
        assert!(residual(&w, &gammas, &sol.dvs) < 1e-5);
        // Fuel = ‖(1,0,0)‖ (Norm) + 1 (vertex gauge) = 2.
        assert_relative_eq!(sol.objective, 2.0, epsilon = 1e-3);
    }

    #[test]
    fn unreachable_w_is_solver_failed() {
        // Γ = [I;0] reaches only coords 0..2; w in coord 4 is unreachable ->
        // the equality is primal-infeasible -> SolverFailed through the wrapper.
        let w = w6([0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
        let gammas = [gamma_top_identity()];
        let err = min_fuel_socp(&w, &gammas, &[norm_gen()]).unwrap_err();
        assert!(matches!(err, PlannerError::SolverFailed(_)));
    }
}
```

- [ ] **Step 5: Declare and re-export the module**

In `src/solver/mod.rs`, add the module declaration and re-export (next to the existing two):

```rust
pub mod extract_qp;
pub mod min_fuel_socp;
pub mod refine_socp;

pub use extract_qp::extract_qp;
pub use min_fuel_socp::{min_fuel_socp, MinFuelSolution};
pub use refine_socp::{refine_socp, RefineSolution};
```

In `src/lib.rs`, extend the solver re-export and the types re-export:

```rust
pub use solver::{extract_qp, min_fuel_socp, refine_socp, MinFuelSolution, RefineSolution};
pub use types::{
    ConicRows, Dual, FuelGenerator, Maneuver, PlannerError, Pseudostate, Solution, SolveParams,
    TimeGrid, M, N,
};
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `cargo test --lib min_fuel_socp 2>&1 | tail -25`
Expected: PASS (all 7 tests). If `single_norm_time_recovers_exact_dv` is off by more than `1e-4`, the budget is not binding here (no budget constraint exists), so this is a clean SOCP — a failure means a cone-assembly bug; re-check the `s = b - A·x` signs.

- [ ] **Step 7: Run fmt + clippy**

Run: `cargo fmt && cargo clippy --all-features -- -D warnings 2>&1 | tail -15`
Expected: clean.

- [ ] **Step 8: Commit**

```bash
git add src/solver/min_fuel_socp.rs src/solver/mod.rs src/lib.rs
git commit -m "feat(solver): add min_fuel_socp — direct full-DOF min-fuel SOCP for robust extraction"
```

---

## Task 3: Rewire Algorithm 3 (`extract`) onto `min_fuel_socp`

Switch `extract` to recover full-DOF maneuvers via `min_fuel_socp`, prune near-zero maneuvers, and update the `solve` call site. `extract_qp` stays (Decision D4) but leaves the `solve` path.

**Files:**
- Modify: `src/algorithm/extract.rs` (rewrite `extract`; update its 2 unit tests)
- Modify: `src/algorithm/mod.rs` (update the `extract::extract(...)` call)

**Interfaces:**
- Consumes: `min_fuel_socp`, `MinFuelSolution` (Task 2); `SublevelSet::fuel_generator` (Task 1).
- Produces: `fn extract<C: CostModel>(cost: &C, grid: &TimeGrid, gammas: &[SMatrix<f64, N, M>], w: &Pseudostate, budget: f64, t_opt: &[usize]) -> Result<ExtractOutcome, PlannerError>` — `ExtractOutcome { maneuvers, total_dv, residual }` unchanged.

- [ ] **Step 1: Update the two existing `extract` unit tests to the new signature**

In `src/algorithm/extract.rs` `mod tests`, the calls `extract(&cost, &grid, &gammas, &w, &q, &lambda, 13.0, &[0])` lose `&q` and `&lambda`. Replace the two test bodies' calls:

```rust
        // extract_recovers_single_maneuver_with_zero_residual:
        let out = extract(&cost, &grid, &gammas, &w, 13.0, &[0]).unwrap();
```

```rust
        // extract_drops_zero_support_times -> rename intent: the t=8 time has
        // Γ=0, so min-fuel assigns it Δv≈0 and pruning drops it.
        let out = extract(&cost, &grid, &gammas, &w, 13.0, &[0, 8]).unwrap();
```

Remove the now-unused `let q = SMatrix::<f64, N, N>::identity();` and `let lambda = ...;` bindings from both tests. Keep all the assertions (they still hold: 1 maneuver, dv ≈ (3,4,12), total_dv ≈ 13, residual < 1e-3, and the second test still ends with 1 maneuver at t=0). Rename the second test to `extract_drops_unused_times` to match the new mechanism, and update its comment.

- [ ] **Step 2: Run to verify failure**

Run: `cargo test --lib algorithm::extract 2>&1 | tail -20`
Expected: FAIL — `extract` still has the old 8-arg signature, arity mismatch.

- [ ] **Step 3: Rewrite `extract`**

Replace the top of `src/algorithm/extract.rs` (imports, consts, and the `extract` fn — keep `ExtractOutcome` as is) with:

```rust
//! Algorithm 3 - Control-Input Extraction: a direct min-fuel SOCP over the
//! converged active set `T^opt` recovers full 3-DOF maneuvers (robust on the
//! degenerate flat contacts where the fixed-support-direction QP under-spans w).

use crate::cost::CostModel;
use crate::solver::min_fuel_socp;
use crate::types::{FuelGenerator, Maneuver, PlannerError, Pseudostate, TimeGrid, M, N};
use nalgebra::{SMatrix, SVector};

/// Maneuvers whose magnitude is below this fraction of the largest recovered
/// maneuver are interior-point dust and are pruned from the reported plan.
/// (Interior-point solutions are not exactly sparse; see Decision D5.)
const PRUNE_REL: f64 = 1e-3;

/// Result of Algorithm 3.
#[derive(Debug, Clone)]
pub(super) struct ExtractOutcome {
    pub maneuvers: Vec<Maneuver>,
    pub total_dv: f64,
    pub residual: f64,
}

/// Algorithm 3 — recover maneuvers over `T^opt` by direct min-fuel SOCP.
///
/// `budget` is the dual optimum `c*` from refinement; it is used only as a
/// self-consistency sanity reference (the SOCP objective should match it to
/// solver tolerance), not as a constraint.
pub(super) fn extract<C: CostModel>(
    cost: &C,
    grid: &TimeGrid,
    gammas: &[SMatrix<f64, N, M>],
    w: &Pseudostate,
    budget: f64,
    t_opt: &[usize],
) -> Result<ExtractOutcome, PlannerError> {
    // Build the per-candidate-time dynamics and fuel generators over T^opt.
    let gammas_t: Vec<SMatrix<f64, N, M>> = t_opt.iter().map(|&k| gammas[k]).collect();
    let generators: Vec<FuelGenerator> = t_opt
        .iter()
        .map(|&k| cost.at(grid.time(k)).fuel_generator())
        .collect();

    let sol = min_fuel_socp(w, &gammas_t, &generators)?;
    debug_assert!(
        budget <= 0.0 || (sol.objective - budget).abs() / budget < 5e-2,
        "min-fuel objective {} disagrees with dual budget {budget}",
        sol.objective
    );

    // Residual of the FULL (unpruned) min-fuel solution: this is the ~0 we report.
    let mut w_acc_full = SVector::<f64, N>::zeros();
    for (idx, &k) in t_opt.iter().enumerate() {
        w_acc_full += gammas[k] * sol.dvs[idx];
    }
    let residual = (w - w_acc_full).norm() / w.norm();

    // Prune interior-point dust: keep maneuvers >= PRUNE_REL of the largest.
    let max_dv = sol
        .dvs
        .iter()
        .map(|dv| dv.norm())
        .fold(0.0_f64, f64::max);
    let keep = PRUNE_REL * max_dv;
    let mut maneuvers = Vec::new();
    let mut total_dv = 0.0;
    for (idx, &k) in t_opt.iter().enumerate() {
        let dv = sol.dvs[idx];
        if dv.norm() <= keep {
            continue;
        }
        total_dv += dv.norm();
        maneuvers.push(Maneuver { t: grid.time(k), dv });
    }

    Ok(ExtractOutcome {
        maneuvers,
        total_dv,
        residual,
    })
}
```

Keep the `#[cfg(test)] mod tests` below (with the Step-1 edits). Note the unused-import lints: `SVector` is still used (`w_acc_full`); `N`, `M` still used.

- [ ] **Step 4: Update the `solve` call site**

In `src/algorithm/mod.rs`, the `extract::extract(...)` call drops `&params.q` and `&refined.lambda`:

```rust
    // --- Algorithm 3: control-input extraction. ---
    let extracted = extract::extract(
        cost,
        &grid,
        &gammas,
        &w,
        refined.objective,
        &refined.t_opt,
    )?;
```

(`params.q` is no longer read here. `SolveParams.q` remains a public field — leave it; add a one-line doc note in `src/types.rs` on the `q` field: `/// (Reserved; unused by the Phase-5b min-fuel extractor.)` appended to its existing doc comment.)

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cargo test --lib algorithm 2>&1 | tail -25`
Expected: PASS — both `extract` unit tests and all `refine`/`init` tests.

- [ ] **Step 6: Run the full integration + gate**

Run: `cargo test --all-features 2>&1 | tail -30`
Expected: `tests/algorithm.rs` still green (`solve_converges_on_reachable_synthetic_problem` now uses min-fuel — its `residual < 1e-2` assert should now be comfortably met, likely < 1e-5).
Run: `cargo fmt && cargo clippy --all-features -- -D warnings 2>&1 | tail -15`
Expected: clean. (If `extract_qp`'s import in `extract.rs` is gone, confirm `extract_qp` itself is still referenced by `lib.rs`/`solver/mod.rs` re-exports + its own tests so no `dead_code` fires.)

- [ ] **Step 7: Commit**

```bash
git add src/algorithm/extract.rs src/algorithm/mod.rs src/types.rs
git commit -m "feat(algorithm): extract via min_fuel_socp (full-DOF, robust on degenerate contacts)"
```

---

## Task 4: Worked-example acceptance (residual collapses)

Tighten the §7 self-consistency test now that extraction recovers `w`, and update the `mdot` report. **Characterization-first** (Decision D5): observe before locking.

**Files:**
- Modify: `examples/mdot.rs` (the self-check at the bottom; reporting unchanged otherwise)
- Modify: `tests/worked_example.rs` (tighten `worked_example_is_self_consistent`)

**Interfaces:**
- Consumes: the public `solve` (unchanged signature) and the now-robust extraction.

- [ ] **Step 1: Characterize — run the example and read the numbers**

Run: `cargo run --example mdot 2>&1 | tail -25`
Record: the printed `residual w_err/w`, `total dv (computed)`, the maneuver count, and the per-maneuver `|dv|` column. Expectations to confirm: residual now ≈ 1e-6–1e-4 (was ~4e-3); total dv ≈ 80.9 mm/s; a *small* maneuver count (single digits) after `PRUNE_REL` pruning. If the count is large (say > 15), inspect the `|dv|` spread and, if there is a clear gap, raise `PRUNE_REL` (e.g. to `1e-2`) in `src/algorithm/extract.rs` and re-run; commit that tuning as part of Task 3's intent (amend or a follow-up `fix` commit). Do NOT loosen the residual assertion to hide a real miss.

- [ ] **Step 2: Write the tightened test**

In `tests/worked_example.rs`, replace the final three asserts of `worked_example_is_self_consistent` (the `total_dv`, `residual`, and `lambda` block, currently asserting `residual < 5e-2`) with:

```rust
    // Phase 5b: the direct min-fuel SOCP now recovers w to ~0 residual with a
    // small maneuver set (the fixed-support QP previously left ~0.4% over ~9
    // maneuvers). Bands set from the characterized run (Step 1).
    assert!(
        sol.total_dv > 0.078 && sol.total_dv < 0.083,
        "total_dv = {} (expected ~80.9 mm/s)",
        sol.total_dv
    );
    assert!(
        sol.residual < 1e-3,
        "residual = {:.3e} (Phase 5b target: << 0.1%)",
        sol.residual
    );
    assert!(
        sol.maneuvers.len() <= 12,
        "expected a small maneuver set, got {}",
        sol.maneuvers.len()
    );
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
```

(Set the `<= 12` bound and the `< 1e-3` residual from the Step-1 characterization — both should hold with margin. Keep `paper_table_iv_does_not_reconstruct` untouched.)

- [ ] **Step 3: Update the `mdot` self-check**

In `examples/mdot.rs`, the bottom self-check currently only asserts iterations and `max_g`. Add a residual assertion so the example is also a guard:

```rust
    assert!(
        sol.residual < 1e-3,
        "Phase 5b: extraction residual {:.3e} should be << 0.1%",
        sol.residual
    );
```

(Place it just after the existing `max_g` assert. The narrative comment block at the top of the file is still accurate — the dynamics are FD-verified and the paper figures remain non-reproducible — but append one sentence: `// Phase 5b: the min-fuel SOCP extractor now reconstructs w to ~0 residual.`)

- [ ] **Step 4: Run**

Run: `cargo test --test worked_example 2>&1 | tail -20`
Expected: PASS (`worked_example_is_self_consistent` + `paper_table_iv_does_not_reconstruct`).
Run: `cargo run --example mdot 2>&1 | tail -20`
Expected: residual << 1e-3, small plan; the assertion passes.

- [ ] **Step 5: Commit**

```bash
git add tests/worked_example.rs examples/mdot.rs
git commit -m "test(worked-example): assert min-fuel extraction recovers w (<0.1% residual, small plan)"
```

---

## Task 5: Hunter L2 cross-check acceptance

The second worked example (Hunter & D'Amico 2025) is pure L2 — the cleanest end-to-end check of the new extractor. Target: residual < 0.01%.

**Files:**
- Modify: `tests/worked_example.rs` (add the Hunter test)

**Interfaces:**
- Consumes: `solve`, `refine_socp`, `AbsoluteOrbit::new`, `J2Roe::new`, `Piecewise::new` (all existing public).

- [ ] **Step 1: Characterize the Hunter case**

Add a temporary print harness (or a `#[test]` with `-- --nocapture`) building the Hunter chief and printing `sol.residual`, `sol.total_dv`, `sol.iterations`, and the all-times dual. The chief (spec §7 second example): `a = 25000 km`, `e_x = -0.658`, `e_y = -0.239` ⇒ `e = √(e_x²+e_y²) ≈ 0.7001`, `ω = atan2(e_y, e_x)` (Rust `atan2` returns ≈ −2.7932 rad = −160°, i.e. **200° mod 360**), `M = u₀ - ω` with `u₀ = 65°` evaluates to ≈ **+3.9277 rad (= 225° = −135° mod 360)**; `i = 51°`, `Ω = 30°`. Don't worry about the angle wrapping: `AbsoluteOrbit::new` takes `ω`/`M` in radians and the Kepler solver normalizes via `wrap_to_pi`, so the raw `atan2`/subtraction outputs produce the correct orbit. (These angles only ever enter through `cos`/`sin`.) Window `39000 s`, `10 s` grid (3901 times), pure L2 via `Piecewise::new(1e18)`. `w = [0.66,-1.52,-0.38,-1.44,0.29,-0.91] / a_c`. Confirm residual < 1e-4 and total_dv ≈ 2.48e-4 m/s (our optimum; the paper's 2.294e-4 bound is NOT our target — see Global Constraints). Lock bands from the observed numbers.

- [ ] **Step 2: Write the test**

Append to `tests/worked_example.rs`:

```rust
#[test]
fn hunter_l2_cross_check_recovers_w() {
    // Hunter & D'Amico 2025 "Sequential Formulation Validation": identical J2 ROE
    // dynamics, pure L2 cost. The dual lower bound is correct (~2.48e-4 m/s in our
    // FD-verified dynamics; the paper's 2.294e-4 is not reproducible — opposite-
    // sign discrepancy, a paper inconsistency). What we assert is that the Phase-5b
    // min-fuel extractor recovers w to <0.01% residual at the self-consistent dual.
    let e_x: f64 = -0.658;
    let e_y: f64 = -0.239;
    let e = (e_x * e_x + e_y * e_y).sqrt();
    let argp = e_y.atan2(e_x); // atan2 -> -2.7932 rad (-160 deg = 200 deg mod 360)
    let u0 = 65.0_f64.to_radians();
    let mean_anom = u0 - argp; // u0 = argp + M -> +3.9277 rad (= -135 deg mod 360); Kepler wraps it
    let chief = AbsoluteOrbit::new(
        A_C,
        e,
        51.0_f64.to_radians(),
        30.0_f64.to_radians(),
        argp,
        mean_anom,
    );
    let dynamics = J2Roe::new(chief, 0.0, 39_000.0);
    let cost = Piecewise::new(1.0e18); // pure Norm2 (no perigee window ever active)
    let w = Pseudostate::from_row_slice(&[0.66, -1.52, -0.38, -1.44, 0.29, -0.91]) / A_C;
    let grid = TimeGrid::uniform(0.0, 39_000.0, 10.0);
    assert_eq!(grid.len(), 3901);

    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).expect("should solve");

    // Self-consistency: refinement objective equals the exact all-times dual.
    let rows: Vec<_> = grid
        .times()
        .map(|t| cost.at(t).cone_constraints(&dynamics.gamma(t)))
        .collect();
    let exact_dual = refine_socp(&w, &rows).expect("exact SOCP").objective;
    assert!(
        (sol.lambda.dot(&w) - exact_dual).abs() / exact_dual < 1e-2,
        "refine dual {} vs exact {exact_dual}",
        sol.lambda.dot(&w)
    );

    // Phase 5b acceptance: extraction reconstructs w to < 0.01% residual.
    assert!(sol.residual < 1e-4, "residual = {:.3e}", sol.residual);
    // Our FD-verified optimum (~2.48e-4 m/s), NOT the paper's 2.294e-4 bound.
    assert!(
        (2.0e-4..=3.0e-4).contains(&sol.total_dv),
        "total_dv = {:.4e} m/s",
        sol.total_dv
    );
    assert!((1..=50).contains(&sol.iterations));
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
}
```

(Tighten the `total_dv` band to the Step-1 observation if it lands cleanly near 2.48e-4.)

- [ ] **Step 3: Run**

Run: `cargo test --test worked_example hunter 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/worked_example.rs
git commit -m "test(worked-example): Hunter L2 cross-check — min-fuel recovers w to <0.01%"
```

---

## Task 6: `active_set_trace` + real-`J2Roe` drop-then-readd refine test

Close the Phase-4 review deferral: expose the per-iteration active set and prove the refinement loop body (drop + re-add) is exercised on the real ill-conditioned dynamics.

**Files:**
- Modify: `src/algorithm/refine.rs` (add the field, populate it)
- Modify: `tests/algorithm.rs` (add the real-`J2Roe` test)

**Interfaces:**
- Produces: `RefineOutcome.active_set_trace: Vec<Vec<usize>>` (the candidate set `T^est` at the start of each iteration's solve, sorted). `#[allow(dead_code)]` (read only by tests, like `max_g_trace`).

- [ ] **Step 1: Write the failing test**

In `tests/algorithm.rs`, add (it will not compile until the field exists, but the algorithmic assertion is the real check — see Step 3 for the characterization):

```rust
use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use std::f64::consts::TAU;

#[test]
fn refine_on_real_j2roe_runs_multiple_iterations() {
    // The Phase-4 well-conditioned synthetic converges too fast to exercise the
    // drop/add loop body. The real worked-example J2Roe Γ is ill-conditioned
    // (δλ row ~1e3 vs others ~1e-4), so refinement takes several iterations on
    // the degenerate e=0.7 contact. This is a public-API guard that the loop
    // body actually runs on real dynamics.
    const A_C: f64 = 25_000e3;
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let dynamics = J2Roe::new(chief, 0.0, 117_990.0);
    let cost = Piecewise::new(TAU / chief.mean_motion());
    let w = SVector::<f64, N>::from_row_slice(&[50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]) / A_C;
    let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0);

    let sol = solve(&dynamics, &cost, w, grid, &SolveParams::default()).expect("should solve");
    assert!(
        sol.iterations >= 2,
        "real J2Roe refinement should take >= 2 iterations, got {}",
        sol.iterations
    );
}
```

(Note: the *drop-then-readd* observation needs the internal `active_set_trace`, which is `pub(super)` and not visible from `tests/`. The public guard above asserts the loop runs more than once on real dynamics. The drop-then-readd assertion lives in a `refine.rs` unit test — Step 4.)

- [ ] **Step 2: Run to verify the public guard fails or passes**

Run: `cargo test --test algorithm refine_on_real_j2roe 2>&1 | tail -15`
Expected: PASS if real refinement already takes ≥ 2 iters (likely). If it converges in 1, this guard is still meaningful as documentation; lower the bound to `>= 1` and rely on the unit test (Step 4) for the loop-body proof. Record the observed iteration count.

- [ ] **Step 3: Add the `active_set_trace` field**

In `src/algorithm/refine.rs`, extend `RefineOutcome`:

```rust
    /// `max_t g` after each solve — non-increasing; read only by tests.
    #[allow(dead_code)]
    pub max_g_trace: Vec<f64>,
    /// The candidate set `T^est` (sorted grid indices) at the start of each
    /// iteration's solve — read only by tests (proves the drop/add loop runs).
    #[allow(dead_code)]
    pub active_set_trace: Vec<Vec<usize>>,
```

Populate it inside `refine`: declare `let mut active_set_trace: Vec<Vec<usize>> = Vec::new();` next to `max_g_trace`, push a snapshot right after the candidate set is fixed for this iteration (just before assembling `rows`):

```rust
        // Solve eq. 40 over the current candidate set T^est.
        let mut snapshot = t_est.clone();
        snapshot.sort_unstable();
        active_set_trace.push(snapshot);
        let rows: Vec<ConicRows> = t_est
            .iter()
            .map(|&k| cost.at(grid.time(k)).cone_constraints(&gammas[k]))
            .collect();
```

Add `active_set_trace` to BOTH `RefineOutcome { … }` constructions (the converged return and — there is only the one converged `Ok(...)`; the `NotConverged` path returns `Err`, so only the success constructor needs the field).

- [ ] **Step 4: Add the drop-then-readd unit test**

In `src/algorithm/refine.rs` `mod tests`, add a test on the **real** `J2Roe` worked example (the `tests` module currently uses only the `SpinDyn` mock — add the real dynamics import locally in the test):

```rust
    #[test]
    fn real_j2roe_refine_drops_then_readds_a_time() {
        use crate::dynamics::{AbsoluteOrbit, J2Roe};
        use std::f64::consts::TAU;
        const A_C: f64 = 25_000e3;
        let chief = AbsoluteOrbit::new(
            A_C,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        let dynamics = J2Roe::new(chief, 0.0, 117_990.0);
        let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0);
        let gammas: Vec<SMatrix<f64, N, M>> = grid.times().map(|t| dynamics.gamma(t)).collect();
        let cost = Piecewise::new(TAU / chief.mean_motion());
        let w = SVector::<f64, N>::from_row_slice(&[50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]) / A_C;
        let params = SolveParams::default();

        // Seed with the coarse init Algorithm 1 would pick, then let refine run.
        let t_est = crate::algorithm::init::initialize(&cost, &grid, &gammas, &w, &params);
        let out = refine(&cost, &grid, &gammas, &w, &params, t_est, 50).unwrap();

        // Characterization aid (run with `-- --nocapture`):
        //   for (i, s) in out.active_set_trace.iter().enumerate() { eprintln!("iter {i}: {s:?}"); }
        assert!(out.iterations >= 2, "iterations = {}", out.iterations);

        // A time dropped in one iteration and present again later == drop-then-readd.
        let trace = &out.active_set_trace;
        let mut readd = false;
        for (i, set_i) in trace.iter().enumerate() {
            for &k in set_i {
                let dropped_later = trace[i + 1..].iter().any(|s| !s.contains(&k));
                let present_after_drop = trace[i + 1..]
                    .iter()
                    .skip_while(|s| s.contains(&k))
                    .any(|s| s.contains(&k));
                if dropped_later && present_after_drop {
                    readd = true;
                }
            }
        }
        assert!(
            readd,
            "expected some time to be dropped then re-added across iterations; trace = {trace:?}"
        );
    }
```

If `init::initialize` is not reachable as `crate::algorithm::init::initialize` from the `refine` test module, call the local `super`-visible path instead (the modules are siblings under `algorithm`; use `super::super::init::initialize`). Confirm visibility by compiling.

> **Characterization-first caveat (Decision D5):** the drop-then-readd assertion is the *goal*, not a guarantee of the default seed. Run with `-- --nocapture` and read the trace first. If the default `n_init`/grid does not exhibit a re-add, adjust the seed (e.g. start from a deliberately suboptimal `t_est` like `vec![0, grid.len()/2, grid.len()-1]`) until the loop demonstrably drops and re-adds, then lock that seed. Do not weaken the assertion to pass — the point is to prove the loop body integrates.

- [ ] **Step 5: Run**

Run: `cargo test --lib algorithm::refine 2>&1 | tail -20`
Run: `cargo test --test algorithm 2>&1 | tail -20`
Expected: PASS (after seed tuning per the caveat).

- [ ] **Step 6: clippy (the new field must not warn)**

Run: `cargo clippy --all-features -- -D warnings 2>&1 | tail -15`
Expected: clean (the `#[allow(dead_code)]` covers `active_set_trace`).

- [ ] **Step 7: Commit**

```bash
git add src/algorithm/refine.rs tests/algorithm.rs
git commit -m "test(refine): expose active_set_trace; prove drop-then-readd on real J2Roe"
```

---

## Task 7: Spec update, full gate, PR

**Files:**
- Modify: `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`

- [ ] **Step 1: Mark Phase 5b done in the spec**

In §6, change the `### Phase 5b` heading from `(fresh-session handoff)` to `✅ Done` and add a short closing paragraph after the existing handoff text:

```markdown
**✅ Done (2026-06-18).** Algorithm 3 now extracts via a direct min-fuel SOCP
(`src/solver/min_fuel_socp.rs`): `min Σⱼ f_{tⱼ}(Δvⱼ) s.t. Σⱼ Γ(tⱼ)Δvⱼ = w` over
`T^opt`, recovering full 3-DOF maneuvers charged by the true per-time cost (L2 SOC
for Norm2 times; a `V_vertex` nonnegative-combination LP for FaceMax times, via the
new `SublevelSet::fuel_generator`). The worked example now reconstructs `w` to
< 0.1 % residual with a small maneuver set, and the Hunter L2 case to < 0.01 %
(`tests/worked_example.rs`), both self-consistent with the exact all-times dual.
The faithful Algorithm-3-as-printed QP (`extract_qp`) is retained as a primitive but
is off the `solve` path. The Phase-4 deferral is closed: `RefineOutcome.active_set_trace`
plus a real-`J2Roe` drop-then-readd refine test.
```

Update the top-of-file `Status:` line: `Phases 0–5b complete. Resume at Phase 6 (Monte Carlo).`

- [ ] **Step 2: Full gate**

Run: `cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features 2>&1 | tail -30`
Expected: all green. Confirm the dynamics guards `tests/fd_stm.rs` and `tests/fd_b_matrix.rs` still pass (they must — no dynamics changed).

- [ ] **Step 3: Commit + push + PR**

```bash
git add docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md
git commit -m "docs(spec): Phase 5b done — robust min-fuel extraction"
git push -u origin phase5b-robust-extraction
gh pr create --title "Phase 5b — robust Algorithm-3 extraction (direct min-fuel SOCP)" \
  --body "Replaces the fixed-support magnitude QP with a direct min-fuel SOCP over T^opt; worked example <0.1% residual, Hunter L2 <0.01%. Dynamics untouched (FD guards green). Closes the Phase-4 active-set-trace deferral."
```

---

## Self-Review (run after implementing, against spec §6 Phase 5b)

1. **Spec coverage:** Candidate fix #1 (direct min-fuel SOCP) → Tasks 1–3. Acceptance "worked example < 0.1%" → Task 4; "Hunter L2 < 0.01%" → Task 5. "Also deferred" (`active_set_trace` + drop-then-readd) → Task 6. "Keep the dual self-consistency assertion" → Tasks 4/5 retain the `refine_socp` all-times comparison. "Don't touch dynamics" → Global Constraints + Task 7 Step 2 re-runs the FD guards. ✓
2. **Placeholder scan:** none — every code/test block is complete; thresholds (`PRUNE_REL`, residual/count bands) are explicitly characterization-first with a stated default and a tuning procedure.
3. **Type consistency:** `FuelGenerator` (Task 1) ↔ `min_fuel_socp` signature (Task 2) ↔ `extract` generators (Task 3); `MinFuelSolution.{dvs,objective}` consumed in Task 3; `extract` arity (6 args) matches the Task-4 call site; `RefineOutcome.active_set_trace: Vec<Vec<usize>>` defined in Task 6 Step 3 and read in Step 4.
