# Phase 3 — Solver Wrappers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wrap the `clarabel` conic solver in two stateless, pure functions — `refine_socp` (Koenig eq. 40 refinement SOCP) and `extract_qp` (Algorithm 3 extraction QP) — that the Phase 4 orchestration will call.

**Architecture:** Both functions translate the paper's math into clarabel's standard cone form `minimize ½xᵀPx + qᵀx s.t. Ax + s = b, s ∈ K`, run `DefaultSolver`, map the terminal `SolverStatus` to `Result`, and return native nalgebra/`Vec` types. They consume *pre-assembled* data (`ConicRows` from the cost layer; `yⱼ` vectors from the caller) so they hold no reference to `Dynamics`/`CostModel` and are unit-testable with hand-built inputs against closed-form optima.

**Tech Stack:** Rust 2021, `nalgebra 0.35` (`SVector`/`SMatrix`), `clarabel 0.11.1` (`CscMatrix`, `DefaultSolver`, `SupportedConeT`), `thiserror 2.0`, `approx 0.5` (dev).

## Global Constraints

- **Work on a branch:** create `phase3-solver-wrappers` off `main` before Task 1; open a PR at the end (mirrors Phase 1/2: `phase2-cost-models` → PR #1).
- **CI gate (run after every task, must be green):**
  `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
- **Dimensions are fixed:** `N = 6` (state/dual `λ`), `M = 3` (control). Already in `src/types.rs`.
- **clarabel solves** `minimize ½xᵀPx + qᵀx  s.t.  Ax + s = b,  s ∈ K`. Note `s = b − Ax`. Cones: `ZeroConeT(n)` (`s=0`), `NonnegativeConeT(n)` (`s ≥ 0` elementwise), `SecondOrderConeT(n)` (`s₀ ≥ ‖s₁…ₙ₋₁‖₂`). Imports: `use clarabel::algebra::*;` (CscMatrix) and `use clarabel::solver::*;` (DefaultSolver, DefaultSettings, DefaultSettingsBuilder, SolverStatus, cone constructors).
- **`P` must be passed UPPER-TRIANGULAR.** clarabel reads only the upper triangle of `P` and does **not** symmetrize (`kkt_assembly.rs`: *"user provided P is always triu regardless"*). Build `P` with the strict lower triangle zeroed; `CscMatrix::from(&dense)` drops exact zeros (`csc/core.rs` filters `v != 0`).
- **`CscMatrix::from(&rows)`** accepts `&Vec<Vec<f64>>` (or `&[[f64; n]]`): outer iter = rows, inner = columns; it asserts all rows have equal length, sets `m = rows.len()`, and stores only nonzeros (all-zero rows still count toward `m`). Use it for `A` and (triu) `P`. Use `CscMatrix::<f64>::zeros((m, n))` for the all-zero `P` of the SOCP.
- **Success statuses:** accept `SolverStatus::Solved` **and** `SolverStatus::AlmostSolved` (reduced accuracy — accepting it keeps Phase 4 refinement robust; the `ε = 0.01` bands absorb the lower accuracy). Every other variant maps to `PlannerError::SolverFailed` with a message naming the status.
- **Assert on bands, not bit-equality** (Risk R3): objective/pinned coordinates to `1e-6`; free coordinates `< 1e-4`; tie expected values to `√` expressions, never decimal literals (avoids `clippy::approx_constant` — bit Phase 2 on `1/√2`).

## Design Decisions (locked, with rationale)

1. **`refine_socp` is a pure SOCP wrapper** taking `(w, &[ConicRows])` → `RefineSolution { lambda, objective }`. It does **not** take `Dynamics`/`CostModel`. The caller (Phase 4 `refine.rs`) assembles `Vec<ConicRows>` via `cost.at(t).cone_constraints(&dynamics.gamma(t))` for each `t ∈ T^est`.

2. **"Per-time slack" is the caller's job, not `refine_socp`'s** (resolves the §4.2 / Phase-3 spec wording "return λ, objective, per-time slack"). Rationale: (a) `refine_socp` consumes `ConicRows`, not cost objects, so it *cannot* compute the contact `g_{U(1,t)}(Γᵀ(t)λ)`; (b) Algorithm 2 must scan `g` over the **full grid `T`** (not just `T^est`) to find local maxima and the global max anyway, so computing it inside `refine_socp` would duplicate that scan. The "slack" Algorithm 2 needs is the **contact value** `g` (comparable to `1`), recomputed via `SublevelSet::contact` — **not** clarabel's raw cone slack `s`. This is documented in the `refine_socp` rustdoc so the Phase 4 author is not surprised.

3. **`λ` is a free (sign-unconstrained) variable.** The SOCP must have **no** sign cone on `λ`. (Guard against copy-pasting the QP's `α ≥ 0`.)

4. **Recover `c*` directly:** `objective = w·λ` from the returned primal `λ`, not from `solver.solution.obj_val` (which is the minimized value `−w·λ`, easy to sign-flip).

5. **`extract_qp` symmetrizes `Q` defensively** (`(Q + Qᵀ)/2`) before forming `P = 2·YᵀQY`, so a non-symmetric `Q` cannot silently corrupt the triu packing.

6. **Warm-starting (Risk R8) is OUT of scope.** Both functions are stateless per call. Phase 4 may revisit; Phase 3 keeps the pure contract.

## File Structure

- **Modify `src/solver/mod.rs`** — add two `pub(crate)` helpers shared by both wrappers (`silent_settings`, `check_status`) and re-export the public functions/types.
- **Modify `src/solver/refine_socp.rs`** — `RefineSolution` struct + `refine_socp` + unit tests (currently a one-line doc-comment stub).
- **Modify `src/solver/extract_qp.rs`** — `extract_qp` + unit tests (currently a one-line doc-comment stub).
- **Modify `src/lib.rs`** — re-export `refine_socp`, `extract_qp`, `RefineSolution` at crate root.
- **Create `tests/solver.rs`** — public-API integration test on a realistic mixed (Norm2 + FaceMax) case.

`ConicRows` (already in `src/types.rs`) is consumed verbatim:
- `linear: Vec<(SVector<f64,6>, f64)>` — each `(a, b)` means `aᵀλ ≤ b`. `FaceMax` emits `(Γvₖ, 1.0)` for `k=1..4`.
- `soc: Vec<(SMatrix<f64,3,6>, f64)>` — each `(G, h)` means `‖Gλ‖₂ ≤ h`, where `G = Γᵀ` (3×6). `Norm2` emits `(Γᵀ, 1.0)`.

---

## Task 1: Shared solver helpers (`solver/mod.rs`)

**Files:**
- Modify: `src/solver/mod.rs`
- Test: `src/solver/mod.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Produces:
  - `pub(crate) fn silent_settings() -> clarabel::solver::DefaultSettings<f64>` — default settings with `verbose = false`.
  - `pub(crate) fn check_status(status: clarabel::solver::SolverStatus) -> Result<(), crate::types::PlannerError>` — `Ok` for `Solved`/`AlmostSolved`, `Err(SolverFailed(...))` otherwise.
  - Re-exports: `pub use refine_socp::{refine_socp, RefineSolution};` and `pub use extract_qp::extract_qp;` (the targets land in Tasks 2–3; add the re-exports here once those items exist, or defer the two `pub use` lines to the end of Task 3 — see Step 5).

- [ ] **Step 1: Write the failing test**

In `src/solver/mod.rs`, append:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use clarabel::solver::SolverStatus;

    #[test]
    fn check_status_accepts_solved_and_almost_solved() {
        assert!(check_status(SolverStatus::Solved).is_ok());
        assert!(check_status(SolverStatus::AlmostSolved).is_ok());
    }

    #[test]
    fn check_status_rejects_failures_naming_the_status() {
        for bad in [
            SolverStatus::PrimalInfeasible,
            SolverStatus::DualInfeasible,
            SolverStatus::MaxIterations,
            SolverStatus::MaxTime,
            SolverStatus::NumericalError,
            SolverStatus::InsufficientProgress,
            SolverStatus::Unsolved,
        ] {
            let err = check_status(bad).unwrap_err();
            // The error message must name the underlying clarabel status,
            // so Phase 4 debugging is not blind.
            assert!(format!("{err}").contains(&format!("{bad:?}")));
        }
    }

    #[test]
    fn silent_settings_are_non_verbose() {
        assert!(!silent_settings().verbose);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib solver::tests 2>&1 | tail -20`
Expected: FAIL — `cannot find function check_status`, `cannot find function silent_settings`.

- [ ] **Step 3: Write minimal implementation**

Replace the entire contents of `src/solver/mod.rs` with:

```rust
//! Convex-solver wrappers around `clarabel`: the refinement SOCP (eq. 40) and
//! the extraction QP (Algorithm 3), plus shared settings/status helpers.

pub mod extract_qp;
pub mod refine_socp;

use crate::types::PlannerError;
use clarabel::solver::{DefaultSettings, DefaultSettingsBuilder, SolverStatus};

/// Default clarabel settings with logging suppressed (keeps the test/CI output
/// clean and the per-iteration SOCP solves quiet during Phase 4 refinement).
pub(crate) fn silent_settings() -> DefaultSettings<f64> {
    DefaultSettingsBuilder::default()
        .verbose(false)
        .build()
        .expect("default clarabel settings are always valid")
}

/// Map a clarabel terminal status to a planner result. `Solved` and
/// `AlmostSolved` (reduced accuracy) are accepted; every other status is a
/// failure whose message names the underlying clarabel status.
pub(crate) fn check_status(status: SolverStatus) -> Result<(), PlannerError> {
    match status {
        SolverStatus::Solved | SolverStatus::AlmostSolved => Ok(()),
        other => Err(PlannerError::SolverFailed(format!(
            "clarabel terminated with status {other:?}"
        ))),
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib solver::tests 2>&1 | tail -20`
Expected: PASS (3 tests). The two `pub mod` lines still point at the unchanged stub files, which compile.

- [ ] **Step 5: Commit**

```bash
git add src/solver/mod.rs
git commit -m "feat(solver): shared silent_settings + check_status helpers"
```

> The `pub use refine_socp::{...}` / `pub use extract_qp::{...}` re-exports are added at the end of Task 3 (Step they reference), once those items exist, to keep this task compiling on its own.

---

## Task 2: `refine_socp` — the eq. 40 refinement SOCP (`solver/refine_socp.rs`)

**Files:**
- Modify: `src/solver/refine_socp.rs`
- Test: `src/solver/refine_socp.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `silent_settings`, `check_status` (Task 1); `ConicRows`, `Dual`, `Pseudostate`, `PlannerError`, `M`, `N` (`crate::types`).
- Produces:
  - `pub struct RefineSolution { pub lambda: Dual, pub objective: f64 }`
  - `pub fn refine_socp(w: &Pseudostate, rows: &[ConicRows]) -> Result<RefineSolution, PlannerError>`

**Encoding (verified against clarabel's shipped `example_socp.rs`):**
`x = λ ∈ ℝ⁶`, `P = 0₆ₓ₆`, `q = −w` (turns `max wᵀλ` into `min −wᵀλ`). Rows assembled **linear-first, then SOC**:
- **Linear** `aᵀλ ≤ b`: one `A`-row `aᵀ`, rhs `b`; all linear rows grouped into one `NonnegativeConeT(n_linear)` (`s = b − aᵀλ ≥ 0`).
- **SOC** `‖Gλ‖₂ ≤ h` (`G` is 3×6): a 4×6 `A`-block `[ 0ᵀ ; −G ]` (first row zeros), rhs `[ h, 0, 0, 0 ]`, one `SecondOrderConeT(4)` (so `s = (h, Gλ)`, and `s₀ ≥ ‖s₁…₃‖` ⟺ `h ≥ ‖Gλ‖`). The `A` and cone vectors are built in the **same** iteration order so blocks align.

- [ ] **Step 1: Write the failing test**

Replace the contents of `src/solver/refine_socp.rs` with the doc comment kept, plus this test module (implementation added in Step 3):

```rust
//! Builds and solves eq. 40 over a candidate-time set (linear + SOC cones from
//! each time's `cone_constraints`), maximize -> minimize.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cost::{FaceMax, Norm2, SublevelSet};
    use crate::types::{ConicRows, M, N};
    use approx::assert_relative_eq;
    use nalgebra::{SMatrix, SVector};

    // Gamma (6x3) with top 3x3 = I, bottom 3x3 = 0. Then Gamma^T = [I_3 | 0],
    // so Gamma^T lambda = (lambda_1, lambda_2, lambda_3); Gamma v = [v; 0].
    fn gamma_top_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(i, i)] = 1.0;
        }
        g
    }

    // Gamma (6x3) with bottom 3x3 = I, top 3x3 = 0. Then Gamma^T = [0 | I_3],
    // so Gamma^T lambda = (lambda_4, lambda_5, lambda_6).
    fn gamma_bottom_identity() -> SMatrix<f64, N, M> {
        let mut g = SMatrix::<f64, N, M>::zeros();
        for i in 0..M {
            g[(M + i, i)] = 1.0;
        }
        g
    }

    fn w6(v: [f64; N]) -> SVector<f64, N> {
        SVector::<f64, N>::from_row_slice(&v)
    }

    // max_t g_{U(1,t)}(Gamma^T(t) lambda) recomputed straight from the rows,
    // normalized by each bound: dual feasibility requires this be <= 1 (+tol).
    fn max_contact(rows: &[ConicRows], lam: &SVector<f64, N>) -> f64 {
        let mut g = f64::NEG_INFINITY;
        for cr in rows {
            for (a, b) in &cr.linear {
                g = g.max(a.dot(lam) / b);
            }
            for (gmat, h) in &cr.soc {
                g = g.max((gmat * lam).norm() / h);
            }
        }
        g
    }

    #[test]
    fn empty_candidate_set_is_invalid_input() {
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let err = refine_socp(&w, &[]).unwrap_err();
        assert!(matches!(err, crate::types::PlannerError::InvalidInput(_)));
    }

    #[test]
    fn s1_single_soc_pins_three_coords() {
        // Norm2 at one time with Gamma^T = [I_3|0]; w hits only those 3 coords.
        // max (3,4,12).(l1,l2,l3) s.t. ||(l1,l2,l3)|| <= 1  ->  ||(3,4,12)|| = 13.
        let rows = vec![Norm2.cone_constraints(&gamma_top_identity())];
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        assert_relative_eq!(sol.objective, 13.0, epsilon = 1e-6);
        assert_relative_eq!(sol.objective, w.dot(&sol.lambda), epsilon = 1e-9);
        assert_relative_eq!(sol.lambda[0], 3.0 / 13.0, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[1], 4.0 / 13.0, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[2], 12.0 / 13.0, epsilon = 1e-6);
        // lambda_4..6 are a free direction (w is zero there): driven ~0.
        for i in 3..N {
            assert!(sol.lambda[i].abs() < 1e-4, "free coord {i} = {}", sol.lambda[i]);
        }
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn s2_face_max_lp_closed_form() {
        // FaceMax (LP path, Risk R6: no published reference) at one time,
        // Gamma=[I_3;0], w=(0,0,1,0,0,0). max l3 s.t. v_k.(l1,l2,l3) <= 1.
        // Binding v3,v4 give b*l3 <= 1 -> l3 = 1/b = sqrt(3); l2 pinned 0.
        let rows = vec![FaceMax.cone_constraints(&gamma_top_identity())];
        let w = w6([0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        let sqrt3 = 3.0_f64.sqrt();
        assert_relative_eq!(sol.objective, sqrt3, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[2], sqrt3, epsilon = 1e-6);
        assert!(sol.lambda[1].abs() < 1e-6);
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn s3_mixed_soc_and_lp_validates_cone_ordering() {
        // The realistic Piecewise case: one FaceMax time (4 linear rows on
        // l1..3) + one Norm2 time (SOC on l4..6). w=(0,0,1,0,0,1).
        // Separable: l3 = sqrt(3) (face-max) and l6 = 1 (||(l4,l5,l6)|| <= 1).
        let rows = vec![
            FaceMax.cone_constraints(&gamma_top_identity()),
            Norm2.cone_constraints(&gamma_bottom_identity()),
        ];
        let w = w6([0.0, 0.0, 1.0, 0.0, 0.0, 1.0]);
        let sol = refine_socp(&w, &rows).unwrap();
        let sqrt3 = 3.0_f64.sqrt();
        assert_relative_eq!(sol.objective, sqrt3 + 1.0, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[2], sqrt3, epsilon = 1e-6);
        assert_relative_eq!(sol.lambda[5], 1.0, epsilon = 1e-6);
        assert!(sol.lambda[1].abs() < 1e-6);
        assert!(sol.lambda[3].abs() < 1e-4);
        assert!(sol.lambda[4].abs() < 1e-4);
        assert!(max_contact(&rows, &sol.lambda) <= 1.0 + 1e-6);
    }

    #[test]
    fn objective_is_scale_equivariant() {
        // Scaling w by k>0 scales c* by k and leaves the lambda direction fixed.
        let rows = vec![Norm2.cone_constraints(&gamma_top_identity())];
        let w = w6([3.0, 4.0, 12.0, 0.0, 0.0, 0.0]);
        let base = refine_socp(&w, &rows).unwrap();
        let scaled = refine_socp(&(w * 2.5), &rows).unwrap();
        assert_relative_eq!(scaled.objective, 2.5 * base.objective, epsilon = 1e-6);
        assert_relative_eq!(scaled.lambda[0], base.lambda[0], epsilon = 1e-6);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib solver::refine_socp 2>&1 | tail -20`
Expected: FAIL — `cannot find function refine_socp`, `cannot find type RefineSolution`.

- [ ] **Step 3: Write minimal implementation**

Insert, between the `//!` doc comment and the `#[cfg(test)]` module, the implementation:

```rust
use crate::solver::{check_status, silent_settings};
use crate::types::{ConicRows, Dual, PlannerError, Pseudostate, M, N};
use clarabel::algebra::CscMatrix;
use clarabel::solver::{DefaultSolver, NonnegativeConeT, SecondOrderConeT, SupportedConeT};

/// The eq. 40 optimum over a candidate-time set: the dual `lambda` and the
/// optimal objective `c* = lambda^T w` (the minimum fuel cost).
///
/// The per-candidate-time contact value `g_{U(1,t)}(Gamma^T(t) lambda)` that
/// Algorithm 2 thresholds against is **not** returned here: it is recomputed by
/// the caller via [`crate::SublevelSet::contact`] (the caller scans `g` over the
/// full grid `T` anyway). Do not confuse it with clarabel's raw cone slack.
#[derive(Debug, Clone)]
pub struct RefineSolution {
    /// Optimal dual `lambda*` in R^6 (outward reachable-set normal).
    pub lambda: Dual,
    /// Optimal objective `c* = lambda*^T w` (>= 0).
    pub objective: f64,
}

/// Solve eq. 40 — `maximize lambda^T w s.t. g_{U(1,t)}(Gamma^T(t) lambda) <= 1`
/// for every candidate time — over `rows`, one [`ConicRows`] per candidate time.
///
/// `lambda` is a free (sign-unconstrained) variable: there is no cone on it.
/// Maps `maximize` to clarabel's `minimize` via `q = -w`; recovers
/// `c* = w . lambda` directly from the primal solution.
pub fn refine_socp(w: &Pseudostate, rows: &[ConicRows]) -> Result<RefineSolution, PlannerError> {
    let n_linear: usize = rows.iter().map(|r| r.linear.len()).sum();
    let n_soc: usize = rows.iter().map(|r| r.soc.len()).sum();
    let total_rows = n_linear + (M + 1) * n_soc;
    if total_rows == 0 {
        return Err(PlannerError::InvalidInput(
            "refine_socp: empty candidate-time set (objective is unbounded)".into(),
        ));
    }

    // Assemble A (total_rows x N) and b in lockstep with the cone vector:
    // all linear rows first (one NonnegativeCone), then one SOC block per SOC.
    let mut a_dense: Vec<Vec<f64>> = Vec::with_capacity(total_rows);
    let mut b: Vec<f64> = Vec::with_capacity(total_rows);
    let mut cones: Vec<SupportedConeT<f64>> = Vec::new();

    if n_linear > 0 {
        cones.push(NonnegativeConeT(n_linear));
        for cr in rows {
            for (a, bi) in &cr.linear {
                a_dense.push(a.iter().copied().collect());
                b.push(*bi);
            }
        }
    }
    for cr in rows {
        for (g, h) in &cr.soc {
            // Scalar-bound row: A all-zero, b = h  ->  s_0 = h.
            a_dense.push(vec![0.0; N]);
            b.push(*h);
            // Vector rows: A = -G, b = 0  ->  s_{1..M} = G lambda.
            for i in 0..M {
                a_dense.push((0..N).map(|c| -g[(i, c)]).collect());
                b.push(0.0);
            }
            cones.push(SecondOrderConeT(M + 1));
        }
    }

    let a_csc = CscMatrix::from(&a_dense);
    let p_csc = CscMatrix::<f64>::zeros((N, N));
    let q: Vec<f64> = (0..N).map(|i| -w[i]).collect();

    let mut solver = DefaultSolver::new(&p_csc, &q, &a_csc, &b, &cones, silent_settings())
        .map_err(|e| PlannerError::SolverFailed(format!("clarabel setup failed: {e:?}")))?;
    solver.solve();
    check_status(solver.solution.status)?;

    let lambda = Dual::from_iterator(solver.solution.x.iter().copied());
    let objective = w.dot(&lambda);
    Ok(RefineSolution { lambda, objective })
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib solver::refine_socp 2>&1 | tail -20`
Expected: PASS (5 tests: `empty_candidate_set_is_invalid_input`, `s1_single_soc_pins_three_coords`, `s2_face_max_lp_closed_form`, `s3_mixed_soc_and_lp_validates_cone_ordering`, `objective_is_scale_equivariant`).

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -25`
Expected: PASS, no clippy warnings. (If clippy flags `needless_range_loop` on an index loop, rewrite with `enumerate`/iterators rather than `#[allow]`.)

- [ ] **Step 6: Commit**

```bash
git add src/solver/refine_socp.rs
git commit -m "feat(solver): refine_socp — eq.40 SOCP via clarabel (SOC + LP cones)"
```

---

## Task 3: `extract_qp` — the Algorithm 3 extraction QP (`solver/extract_qp.rs`)

**Files:**
- Modify: `src/solver/extract_qp.rs`
- Modify: `src/solver/mod.rs` (add the two `pub use` re-exports deferred from Task 1)
- Modify: `src/lib.rs` (crate-root re-exports)
- Test: `src/solver/extract_qp.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `silent_settings`, `check_status` (Task 1); `Pseudostate`, `PlannerError`, `N` (`crate::types`).
- Produces:
  - `pub fn extract_qp(w: &Pseudostate, ys: &[SVector<f64, N>], q_weight: &SMatrix<f64, N, N>, budget: f64) -> Result<Vec<f64>, PlannerError>` — returns the nonnegative magnitudes `αⱼ`, one per `yⱼ`. The caller forms `Maneuver { t: tⱼ, dv: αⱼ·sⱼ }` (Phase 4), where `yⱼ = Γ(tⱼ)·sⱼ` and `sⱼ = support(...)`.

**Encoding (verified):** minimize `(w − Yα)ᵀQ(w − Yα)` with `Y = [y₁|…|y_K]` (6×K). Expand to `αᵀ(YᵀQY)α − 2(YᵀQw)ᵀα + wᵀQw`, drop the constant. clarabel form: `P = 2·YᵀQY` (passed **upper-triangular**), `q = −2·YᵀQw`. Constraints `α ≥ 0` and `Σα ≤ budget` combine into one `NonnegativeConeT(K+1)` with `A = [ −I_K ; 1ᵀ ]`, `b = [ 0_K ; budget ]`.

- [ ] **Step 1: Write the failing test**

Replace the contents of `src/solver/extract_qp.rs` with the doc comment plus this test module (implementation added in Step 3):

```rust
//! The Algorithm 3 QP: solve for nonnegative magnitudes that minimize the
//! weighted pseudostate residual.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::N;
    use approx::assert_relative_eq;
    use nalgebra::{SMatrix, SVector};

    fn e(i: usize) -> SVector<f64, N> {
        let mut v = SVector::<f64, N>::zeros();
        v[i] = 1.0;
        v
    }
    fn w6(v: [f64; N]) -> SVector<f64, N> {
        SVector::<f64, N>::from_row_slice(&v)
    }
    // Weighted residual (w - Y alpha)^T Q (w - Y alpha).
    fn weighted_obj(
        w: &SVector<f64, N>,
        ys: &[SVector<f64, N>],
        q: &SMatrix<f64, N, N>,
        alpha: &[f64],
    ) -> f64 {
        let mut werr = *w;
        for (a, y) in alpha.iter().zip(ys) {
            werr -= *a * *y;
        }
        (werr.transpose() * q * werr)[(0, 0)]
    }

    #[test]
    fn no_directions_is_invalid_input() {
        let w = w6([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let q = SMatrix::<f64, N, N>::identity();
        let err = extract_qp(&w, &[], &q, 10.0).unwrap_err();
        assert!(matches!(err, crate::types::PlannerError::InvalidInput(_)));
    }

    #[test]
    fn qp_a_interior_optimum_exact_fit() {
        // y1=e1, w=2 e1, budget slack -> alpha=2, residual 0.
        let q = SMatrix::<f64, N, N>::identity();
        let a = extract_qp(&w6([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]), &[e(0)], &q, 10.0).unwrap();
        assert_eq!(a.len(), 1);
        assert_relative_eq!(a[0], 2.0, epsilon = 1e-6);
    }

    #[test]
    fn qp_b_budget_binds() {
        // Same as A but budget=1 -> alpha=1, residual 1.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let a = extract_qp(&w, &[e(0)], &q, 1.0).unwrap();
        assert_relative_eq!(a[0], 1.0, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &[e(0)], &q, &a), 1.0, epsilon = 1e-6);
    }

    #[test]
    fn qp_c_nonneg_binds() {
        // w = -3 e1, only nonneg direction available -> alpha=0, residual 3.
        let q = SMatrix::<f64, N, N>::identity();
        let a = extract_qp(&w6([-3.0, 0.0, 0.0, 0.0, 0.0, 0.0]), &[e(0)], &q, 10.0).unwrap();
        assert!(a[0].abs() < 1e-6);
    }

    #[test]
    fn qp_d_two_orthonormal_budget_binds() {
        // y1=e1,y2=e2, w=(2,3,..), budget=4 (sum unconstrained = 5).
        // Equal-shrink -> alpha=(1.5,2.5), weighted obj 0.5.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let ys = [e(0), e(1)];
        let a = extract_qp(&w, &ys, &q, 4.0).unwrap();
        assert_relative_eq!(a[0], 1.5, epsilon = 1e-6);
        assert_relative_eq!(a[1], 2.5, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &ys, &q, &a), 0.5, epsilon = 1e-6);
    }

    #[test]
    fn qp_e_weighted_q_budget_binds() {
        // Q = diag(1,4,1,1,1,1), y1=e1,y2=e2, w=(2,3,..), budget=4.
        // KKT: 2(2-a1)=8(3-a2) with a1+a2=4 -> alpha=(1.2,2.8), weighted obj 0.8.
        let mut q = SMatrix::<f64, N, N>::identity();
        q[(1, 1)] = 4.0;
        let w = w6([2.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let ys = [e(0), e(1)];
        let a = extract_qp(&w, &ys, &q, 4.0).unwrap();
        assert_relative_eq!(a[0], 1.2, epsilon = 1e-6);
        assert_relative_eq!(a[1], 2.8, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &ys, &q, &a), 0.8, epsilon = 1e-6);
    }

    #[test]
    fn qp_f_non_orthogonal_directions_exercise_off_diagonal_p() {
        // y1=e1, y2=e1+e2 (non-orthogonal -> P has off-diagonal terms).
        // w=(2,3,..): unconstrained LS gives alpha1<0, so nonneg binds ->
        // alpha=(0, 2.5), residual^2 = 0.5. Catches a triu-P packing bug.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 3.0, 0.0, 0.0, 0.0, 0.0]);
        let y2 = e(0) + e(1);
        let ys = [e(0), y2];
        let a = extract_qp(&w, &ys, &q, 10.0).unwrap();
        assert!(a[0].abs() < 1e-6, "alpha1 should hit the nonneg bound: {}", a[0]);
        assert_relative_eq!(a[1], 2.5, epsilon = 1e-6);
        assert_relative_eq!(weighted_obj(&w, &ys, &q, &a), 0.5, epsilon = 1e-6);
    }

    #[test]
    fn qp_residual_unique_when_p_singular() {
        // Duplicate directions -> Y^T Q Y singular -> alpha non-unique, but the
        // residual w - Y*alpha is unique. Assert on the residual, not alpha.
        let q = SMatrix::<f64, N, N>::identity();
        let w = w6([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let ys = [e(0), e(0)];
        let a = extract_qp(&w, &ys, &q, 10.0).unwrap();
        assert_relative_eq!(a[0] + a[1], 2.0, epsilon = 1e-5);
        assert!(weighted_obj(&w, &ys, &q, &a) < 1e-8);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --all-features --lib solver::extract_qp 2>&1 | tail -20`
Expected: FAIL — `cannot find function extract_qp`.

- [ ] **Step 3: Write minimal implementation**

Insert, between the `//!` doc comment and the `#[cfg(test)]` module:

```rust
use crate::solver::{check_status, silent_settings};
use crate::types::{PlannerError, Pseudostate, N};
use clarabel::algebra::CscMatrix;
use clarabel::solver::{DefaultSolver, NonnegativeConeT, SupportedConeT};
use nalgebra::{SMatrix, SVector};

/// Algorithm 3 QP: pick nonnegative magnitudes `alpha_j >= 0` (one per maneuver
/// direction) minimizing the `Q`-weighted residual `(w - sum_j alpha_j y_j)^T Q
/// (w - sum_j alpha_j y_j)`, subject to `sum_j alpha_j <= budget`
/// (`budget = lambda_opt^T w`).
///
/// `ys[j] = Gamma(t_j) . s_j` is the pseudostate contribution of the unit
/// support direction `s_j` at the j-th optimal time; the caller builds the
/// `Maneuver` as `dv = alpha_j . s_j` applied at `t_j`.
pub fn extract_qp(
    w: &Pseudostate,
    ys: &[SVector<f64, N>],
    q_weight: &SMatrix<f64, N, N>,
    budget: f64,
) -> Result<Vec<f64>, PlannerError> {
    let k = ys.len();
    if k == 0 {
        return Err(PlannerError::InvalidInput(
            "extract_qp: no maneuver directions".into(),
        ));
    }
    if budget < 0.0 {
        return Err(PlannerError::InvalidInput(format!(
            "extract_qp: budget must be non-negative, got {budget}"
        )));
    }

    // Symmetrize Q defensively so the triu(P) packing cannot drop an
    // asymmetric part. Q is PD (identity by default), so qsym is PD.
    let qsym = (q_weight + q_weight.transpose()) * 0.5;
    let qy: Vec<SVector<f64, N>> = ys.iter().map(|y| qsym * y).collect();
    let qw = qsym * w;

    // P = 2 Y^T Q Y, emitted UPPER-TRIANGULAR (strict-lower zeroed; CscMatrix
    // drops the zeros). Keep the full diagonal value (do not halve it).
    let mut p_dense: Vec<Vec<f64>> = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in i..k {
            p_dense[i][j] = 2.0 * ys[i].dot(&qy[j]);
        }
    }
    let p_csc = CscMatrix::from(&p_dense);

    // q = -2 Y^T Q w.
    let q: Vec<f64> = ys.iter().map(|y| -2.0 * y.dot(&qw)).collect();

    // A = [ -I_K ; 1^T ]  ((K+1) x K), b = [ 0_K ; budget ], one NonnegativeCone.
    let mut a_dense: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
    for i in 0..k {
        let mut row = vec![0.0; k];
        row[i] = -1.0;
        a_dense.push(row);
    }
    a_dense.push(vec![1.0; k]);
    let a_csc = CscMatrix::from(&a_dense);

    let mut b = vec![0.0; k + 1];
    b[k] = budget;
    let cones: [SupportedConeT<f64>; 1] = [NonnegativeConeT(k + 1)];

    let mut solver = DefaultSolver::new(&p_csc, &q, &a_csc, &b, &cones, silent_settings())
        .map_err(|e| PlannerError::SolverFailed(format!("clarabel setup failed: {e:?}")))?;
    solver.solve();
    check_status(solver.solution.status)?;

    Ok(solver.solution.x.clone())
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --all-features --lib solver::extract_qp 2>&1 | tail -25`
Expected: PASS (8 tests).

- [ ] **Step 5: Add the deferred re-exports**

In `src/solver/mod.rs`, add after the `pub mod` lines (before the `use` block is fine):

```rust
pub use extract_qp::extract_qp;
pub use refine_socp::{refine_socp, RefineSolution};
```

In `src/lib.rs`, add a crate-root re-export line after the existing `pub use` block:

```rust
pub use solver::{extract_qp, refine_socp, RefineSolution};
```

- [ ] **Step 6: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -25`
Expected: PASS, no warnings.

- [ ] **Step 7: Commit**

```bash
git add src/solver/extract_qp.rs src/solver/mod.rs src/lib.rs
git commit -m "feat(solver): extract_qp — Algorithm 3 QP via clarabel; wire re-exports"
```

---

## Task 4: Public-API integration test (`tests/solver.rs`)

**Files:**
- Create: `tests/solver.rs`

**Interfaces:**
- Consumes: `koenig_planner::{refine_socp, extract_qp, RefineSolution}` (crate-root re-exports from Task 3); `koenig_planner::cost::{FaceMax, Norm2, SublevelSet}`; `koenig_planner::types::{ConicRows, M, N}`.

This test exercises the two wrappers through the **public** API on a mixed Norm2 + FaceMax problem, then feeds the SOCP's optimum into the QP — the same hand-off Phase 4 will wire — confirming the budget `λᵀw` is consistent end-to-end.

- [ ] **Step 1: Write the failing test**

Create `tests/solver.rs`:

```rust
//! Public-API integration tests for the Phase 3 solver wrappers.

use approx::assert_relative_eq;
use koenig_planner::cost::{FaceMax, Norm2, SublevelSet};
use koenig_planner::{extract_qp, refine_socp};
use nalgebra::{SMatrix, SVector};

const N: usize = 6;
const M: usize = 3;

fn gamma_top_identity() -> SMatrix<f64, N, M> {
    let mut g = SMatrix::<f64, N, M>::zeros();
    for i in 0..M {
        g[(i, i)] = 1.0;
    }
    g
}
fn gamma_bottom_identity() -> SMatrix<f64, N, M> {
    let mut g = SMatrix::<f64, N, M>::zeros();
    for i in 0..M {
        g[(M + i, i)] = 1.0;
    }
    g
}

#[test]
fn refine_then_extract_hands_off_through_public_api() {
    // Mixed problem: FaceMax on l1..3 + Norm2 (SOC) on l4..6, w=(0,0,1,0,0,1).
    // c* = sqrt(3)+1 (validated in unit tests); use it as the QP budget.
    let rows = vec![
        FaceMax.cone_constraints(&gamma_top_identity()),
        Norm2.cone_constraints(&gamma_bottom_identity()),
    ];
    let w = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]);

    let refined = refine_socp(&w, &rows).unwrap();
    assert_relative_eq!(refined.objective, 3.0_f64.sqrt() + 1.0, epsilon = 1e-6);
    assert!(refined.objective >= 0.0);

    // Hand off to the QP: two directions that exactly reconstruct w
    // (y1 = e3, y2 = e6); budget = c* is slack, so alpha = (1, 1), residual 0.
    let y1 = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
    let y2 = SVector::<f64, N>::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);
    let q = SMatrix::<f64, N, N>::identity();
    let alpha = extract_qp(&w, &[y1, y2], &q, refined.objective).unwrap();
    assert_relative_eq!(alpha[0], 1.0, epsilon = 1e-6);
    assert_relative_eq!(alpha[1], 1.0, epsilon = 1e-6);

    let werr = w - alpha[0] * y1 - alpha[1] * y2;
    assert!(werr.norm() < 1e-6);
}
```

- [ ] **Step 2: Run test to verify it fails (then passes)**

Run: `cargo test --all-features --test solver 2>&1 | tail -20`
Expected: with Tasks 1–3 complete this should **pass** immediately (the implementations already exist). If it fails to compile on a missing re-export, fix the Task-3 Step-5 re-exports. (The test is "failing-first" only in the sense that the file did not exist before.)

- [ ] **Step 3: Run the full gate**

Run: `cargo fmt --all -- --check && cargo clippy --all-features -- -D warnings && cargo test --all-features 2>&1 | tail -30`
Expected: PASS. Phase 3 adds **17 tests** to the current suite: 3 (Task 1) + 5 (Task 2) + 8 (Task 3) + 1 (Task 4). (The Phase-2 suite was 44; if the harness counts integration tests the same way, expect ~61 total — but treat the +17 as the contract, not the absolute number.)

- [ ] **Step 4: Commit**

```bash
git add tests/solver.rs
git commit -m "test(solver): public-API integration — refine_socp -> extract_qp hand-off"
```

- [ ] **Step 5: Push and open the PR**

```bash
git push -u origin phase3-solver-wrappers
gh pr create --fill --title "Phase 3 — solver wrappers (refine_socp + extract_qp)" \
  --body "Implements Phase 3 per docs/superpowers/plans/2026-06-17-koenig-planner-phase3-solver-wrappers.md. Closed-form SOCP/QP unit tests (incl. the Risk-R6 face-max LP path and a non-orthogonal-Y triu-P case) + public-API integration test. CI green."
```

---

## Spec & Risk Coverage (self-review)

- **Phase 3 deliverable "refine_socp: assemble eq. 40 ... map maximize→minimize, return λ, objective, per-time slack"** → Task 2. `maximize→minimize` via `q = −w` (verified). "Per-time slack" resolved (Design Decision 2): returned as the caller-recomputed contact, documented in `RefineSolution` rustdoc; `refine_socp` returns `{lambda, objective}`.
- **Phase 3 deliverable "extract_qp: the Algorithm 3 QP"** → Task 3 (`P = 2YᵀQY` triu, `q = −2YᵀQw`, `α ≥ 0`, `Σα ≤ λᵀw`).
- **Exit criterion "solver tests pass on small hand-checkable problems with closed-form optima"** → Tasks 2–4: S1 (pure SOC), S2 (face-max LP), S3 (mixed cone ordering), QP-A…F + singular-P + integration. All optima independently re-derived and confirmed.
- **Risk R2 (clarabel convention mismatch)** → SOC sign/layout verified against clarabel's `example_socp.rs`; cone ordering assembled in lockstep with `A`/`b`; `λ` explicitly free; `check_status` enumerates the full `SolverStatus` set (Task 1).
- **Risk R3 (tolerances)** → band assertions (`1e-6` pinned, `1e-4` free), `√` expressions not literals.
- **Risk R6 (face-max LP has no published reference)** → S2 is a standalone hand-derived closed-form LP optimum (`c* = √3`), not entangled with the eq.-49 piecewise case.
- **Risk R8 (clarabel ≠ paper solver; warm-start)** → out of scope, stateless contract kept (Design Decision 6).
- **Edge cases** → empty candidate set / no directions → `InvalidInput` (not panic); `w = 0` → `c* = 0`; singular `P` → assert on the unique residual, not non-unique `α`; non-orthogonal `Y` → exercises off-diagonal/triu `P`; negative budget → `InvalidInput`.
- **Type consistency** → `ConicRows.linear: (SVector<f64,6>, f64)` and `.soc: (SMatrix<f64,3,6>, f64)` consumed verbatim (no double-transpose: `Norm2` already emits `G = Γᵀ`). `Dual = SVector<f64,6>`, `Pseudostate = SVector<f64,6>`. `RefineSolution` named identically across Tasks 2/3/4.

## Notes for Phase 4 (not implemented here)

- Compute the per-time contact `g_{U(1,t)}(Γᵀ(t)λ)` via `cost.at(t).contact(dynamics.gamma(t).transpose() * λ)` over the full grid for the drop (`< 1 − ε_remove`) / add (`> 1`) / converge (`max_t g ≤ 1 + ε_cost`) logic.
- Before `extract_qp`, drop optimal times whose support direction `sⱼ` is ~zero (a `yⱼ = 0` column leaves `αⱼ` unconstrained-but-irrelevant).
- Toy-`Γ` band tolerances may need loosening once real Phase-1 `Γ(t)` (ill-conditioned near `e = 0.7` perigee) is wired in.

---

**Plan complete and saved to `docs/superpowers/plans/2026-06-17-koenig-planner-phase3-solver-wrappers.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

**Which approach?**
