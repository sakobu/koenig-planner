# Koenig Planner — Phase 2 (Cost Models) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the two validation cost models and the time-varying selector so the `SublevelSet`/`CostModel` traits are fully realized: `Norm2` (`‖u‖₂`), `FaceMax` (`max(Vᶠᵃᶜᵉu)`, tetrahedral occulter, eq. 47–48), and `Piecewise` (eq. 49 — FaceMax in 2-hr perigee windows, Norm2 elsewhere).

**Architecture:** Three focused files under `src/cost/`, replacing the Phase 0 stubs. Each cost implements all three `SublevelSet` methods — `contact` (`g(y)`), `support` (`s(y)`), and `cone_constraints` (the Table II solver rows: one SOC for Norm2, four linear rows for FaceMax). `Piecewise` holds a `Norm2` and a `FaceMax` plus the orbit period and window half-width, and `CostModel::at(t)` dispatches between them via the eq. 49 perigee-window test. Correctness rests on the contact/support identity (eq. 23 `λᵀs(λ)=g(λ)`), positive homogeneity (eq. 8 / Property 3), known-direction reference values, a `V_vertex`/`V_face` transcription cross-check (`f(vₖ)=1/9`), a dynamics-free `cone_constraints↔contact` consistency check, and the window-boundary logic. Implementing `cone_constraints` now (rather than deferring to Phase 3) completes the trait, leaves **no `unimplemented!` in the cost layer**, and de-risks the solver wiring — it is pure cost math fully determined by each model.

**Tech Stack:** Rust 2021, `nalgebra` 0.35 (`SVector<f64,3>`, `SMatrix<f64,6,3>`/`<f64,3,6>`), `approx` 0.5 (dev, float-tolerant asserts). No new dependencies.

**Verification status:** The **entire Phase 2 implementation below was written into a throwaway `git worktree` off `HEAD` and run through the full gate before this plan was finalized** — `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test` all pass (**18 new cost tests green**: 9 `FaceMax` + 6 `Norm2` + 2 `Piecewise` + 1 carried-over wiring test, on top of all Phase 0/1 tests). All reference numbers were produced by an **independent pure-Python oracle** (no nalgebra), and the cost-matrix transcription (eq. 47–49, Table II, eq. 23) was **independently re-verified character-by-character against `docs/Planner.pdf` by a separate agent** — entry-for-entry match. Two findings from that pass are folded in: (1) positive homogeneity is the paper's **eq. 8 / Property 3**, *not* eq. 23 (eq. 23 is only the support identity) — cited correctly throughout; (2) the paper does **not** pin down whether the eq. 49 window period is the Keplerian `2π/n` or the J₂-perturbed period, so `Piecewise` takes the period as a constructor argument (the worked example passes `2π/n` ≈ 10.93 hr = 39338.8 s, consistent with the paper's rounded "10.92 hr"). Two bugs the scratch run caught and that are pre-fixed here: **clippy::approx_constant** fires on the literal `0.70710678…` (≈ `1/√2`) — the L2 contact test uses the `(3,4,12)→13` vector instead; and the **exact** window boundary (`|t−center| = 3600 s`) is a floating-point knife-edge, so the boundary test probes ±1 s either side instead of asserting the tie.

## Global Constraints

These apply to every task; values copied verbatim from the design spec (§4.1, §5.2–5.3) and confirmed by the PDF re-verification.

- **Crate location & gate:** crate root is the project root. The binding gate after every task (matches CI exactly): `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test`. Compare floats in tests with `approx` (matrices included), never `assert_eq!` on `f64`.
- **Reuse `crate::types`.** Use `M` (= 3) and `N` (= 6) and `ConicRows` from `crate::types`; never hardcode `3`/`6`. `y = Γᵀ(t)λ ∈ ℝᴹ` is the contact/support argument; `gamma_t: &SMatrix<f64, N, M>` is `Γ(t)` (N×M) in `cone_constraints`.
- **`ConicRows` encoding (already defined in `types.rs`):** `linear: Vec<(SVector<f64,N>, f64)>` is `(a, b)` with `aᵀλ ≤ b`; `soc: Vec<(SMatrix<f64,M,N>, f64)>` is `(G, h)` with `‖Gλ‖₂ ≤ h`. Do not change this type. Norm2 → one SOC row `(Γᵀ, 1.0)`; FaceMax → four linear rows `(Γ vₖ, 1.0)`.
- **Cost models (Table II), exact:** Norm2 — `g(y)=‖y‖₂`, `s(y)=y/‖y‖₂`. FaceMax — `g(y)=maxₖ yᵀwₖ`, `s(y)=argmaxₖ wₖ`, over the columns `wₖ` of `W = [0₃ₓ₁ | V_vertex]` (origin column prepended). The origin column makes `g(y) ≥ 0` (i.e. `g(y)=max(0, maxₖ yᵀvₖ)` and `g(0)=0`); it is vacuous as a cone row (`0ᵀλ ≤ 1`) and is omitted there.
- **`V_vertex` / `V_face` (eq. 47–48), verified against `docs/Planner.pdf`:**
  ```
  V_vertex = [[ √(2/3), −√(2/3),  0,      0     ],
              [ 0,       0,       √(2/3), −√(2/3)],
              [−√(1/3), −√(1/3),  √(1/3),  √(1/3)]]   (3×4: columns are the 4 unit thruster directions)
  V_face = (1/3)·[[−√(2/3), 0,      √(1/3)],
                 [ √(2/3), 0,      √(1/3)],
                 [ 0,     −√(2/3), −√(1/3)],
                 [ 0,      √(2/3), −√(1/3)]]            (4×3: cost = max over the 4 face-rows of V_face·u)
  ```
  Only `V_vertex` (via `W`) enters the algorithm (contact/support/cone). `V_face` is the cost *definition* and appears in the implementation **only** as a test-module cross-check.
- **Identities (cite correctly):** eq. 23 is the support identity `λᵀs(λ) = g(λ)`. Positive homogeneity `g(αy)=α·g(y)` for `α ≥ 0` is **eq. 8 / Property 3**, a separate result. Both must hold for both costs.
- **Piecewise (eq. 49):** `f = max(Vᶠᵃᶜᵉu)` on `T₁`, `‖u‖₂` on `T₂`. `T₁ = { t : |t − (k+0.5)·period| < 1 hr, k ∈ ℤ }` (strict `<`, half-width 3600 s, i.e. 2-hr windows). The centers `(k+0.5)·period` coincide with perigee for the worked example because its chief starts at apogee (`M₀ = 180°`) at `t = 0`. `Piecewise::new` takes `period` (the paper leaves Keplerian-vs-perturbed open; the worked example uses `2π/n`).
- **Complete the trait — no stubs left.** After Phase 2, every `SublevelSet`/`CostModel` method has a real body; no `unimplemented!`/`#[allow(unused_variables)]` remains in `src/cost/`.
- **`Piecewise` drops `Default`.** It now carries fields, so the Phase 0 `#[derive(..., Default)]` and the `&Piecewise` unit-struct usage in `cost/mod.rs`'s wiring test must change to `Piecewise::new(period)`. `Norm2`/`FaceMax` stay unit structs (keep `Default`).

## Reference values (independent pure-Python oracle — used by the tests below)

```
√(2/3) = 0.816496580927726        √(1/3) = 0.5773502691896257

V_vertex: all 4 columns unit-norm; every pairwise column dot = −1/3 (regular tetrahedron).
f(vₖ) = max(V_face · vₖ) = 1/9 for every vertex column   (V_vertex/V_face transcription cross-check)

FaceMax  (g = max over W = [0 | V_vertex]; s = argmax column, ties → lowest index):
  y=( 1, 0, 0):     g = √(2/3) = 0.816497    s = v₀ = ( √(2/3), 0, −√(1/3))      (unique)
  y=( 0, 0, 1):     g = √(1/3) = 0.577350    s = v₂ = (0, √(2/3),  √(1/3))       (tie v₂/v₃ → v₂)
  y=( 0.3, 0.4,0.5): g = 0.615273766966       (eq. 23: y·s = g)
  y=(−0.6, 0.2,−0.9): g = 1.009513190827      s = v₁ = (−√(2/3), 0, −√(1/3))
  y = 0:            g = 0, s = 0 (origin column wins)

Norm2:
  y=(−0.6, 0.2,−0.9): s = (−0.5454545454545454, 0.18181818181818182, −0.8181818181818182), ‖s‖ = 1
  y=(3,4,12): g = 13.0

cone_constraints consistency (fixed, dynamics-free):
  Gamma  = [[1,0,0],[0,1,0],[0,0,1],[0.5,0.5,0],[0.2,−0.3,0.7],[−0.4,0.1,0.6]]   (6×3)
  lambda = [0.5,−0.2,0.8,0.1,−0.6,0.3]   →   Gammaᵀ·lambda = (0.31, 0.06, 0.56)
  Norm2:   ‖Gammaᵀ lambda‖₂            = 0.6428841264178172   (one SOC row, G=Gammaᵀ, h=1)
  FaceMax: maxₖ (Gamma vₖ)ᵀ lambda     = 0.37230594560185404  (four linear rows, h=1)

Window geometry (eq. 49) for the Piecewise tests: period = 40000 s, half_width = 3600 s
  → perigee centers at (k+0.5)·40000; first window = open interval (16400, 23600) s.
```

---

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `src/cost/norm2.rs` | `Norm2`: `contact`/`support`/`cone_constraints` (replaces stub) | Task 1 |
| `src/cost/facemax.rs` | `FaceMax` + `vertex_columns()`; `contact`/`support`/`cone_constraints` (replaces stub) | Task 2 |
| `src/cost/piecewise.rs` | `Piecewise{norm2, facemax, period, half_width}` + `new`/`in_perigee_window`; `CostModel::at` (replaces stub) | Task 3 |
| `src/cost/mod.rs` | unchanged traits/re-exports; **wiring test** updated for `Piecewise::new` | Task 3 |

The three modules are already declared and re-exported in `src/cost/mod.rs` (`pub mod {norm2,facemax,piecewise}; pub use ...`), so no module-declaration steps are needed — each task replaces a stub file's body. Task order follows dependencies: `Norm2` and `FaceMax` are independent; `Piecewise` consumes both.

---

## Task 1: `Norm2` — the L2 cost (`‖u‖₂`)

Replace the `unimplemented!` stub in `src/cost/norm2.rs` with the real contact/support/SOC-row implementation. Simplest cost; no dependencies on the other two.

**Files:**
- Modify: `src/cost/norm2.rs` (replace the stub `impl`; add a `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `super::SublevelSet`, `crate::types::{ConicRows, M, N}`, `nalgebra::{SMatrix, SVector}` (all already imported in the stub).
- Produces: `impl SublevelSet for Norm2` with
  - `fn contact(&self, y: SVector<f64, M>) -> f64` = `‖y‖₂`
  - `fn support(&self, y: SVector<f64, M>) -> SVector<f64, M>` = `y/‖y‖₂` (`0` at `y=0`)
  - `fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows` = one SOC row `(Γᵀ, 1.0)`.

- [ ] **Step 1: Append the failing tests to `src/cost/norm2.rs`**

Add this test module at the end of the stub file (the `unimplemented!` impls are still in place, so these compile but panic at runtime):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn contact_is_the_l2_norm() {
        // (3,4,12) has the clean norm 13 (avoids any constant-like literal).
        assert_relative_eq!(
            Norm2.contact(SVector::<f64, M>::new(3.0, 4.0, 12.0)),
            13.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            Norm2.contact(SVector::<f64, M>::new(1.0, 0.0, 0.0)),
            1.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn support_is_the_unit_direction() {
        let y = SVector::<f64, M>::new(-0.6, 0.2, -0.9);
        let s = Norm2.support(y);
        assert_relative_eq!(
            s,
            SVector::<f64, M>::new(
                -0.5454545454545454,
                0.18181818181818182,
                -0.8181818181818182
            ),
            epsilon = 1e-12
        );
        assert_relative_eq!(s.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn support_of_zero_is_zero() {
        assert_relative_eq!(
            Norm2.support(SVector::<f64, M>::zeros()),
            SVector::<f64, M>::zeros(),
            epsilon = 1e-15
        );
    }

    #[test]
    fn contact_support_identity_eq23() {
        // eq. 23: lambda . s(lambda) = g(lambda).
        for y in [
            SVector::<f64, M>::new(0.3, 0.4, 0.5),
            SVector::<f64, M>::new(-0.6, 0.2, -0.9),
        ] {
            assert_relative_eq!(y.dot(&Norm2.support(y)), Norm2.contact(y), epsilon = 1e-12);
        }
    }

    #[test]
    fn positive_homogeneity() {
        // g(a y) = a g(y) for a >= 0 (eq. 8 / Property 3).
        let y = SVector::<f64, M>::new(0.3, 0.4, 0.5);
        assert_relative_eq!(
            Norm2.contact(y * 2.5),
            2.5 * Norm2.contact(y),
            epsilon = 1e-12
        );
    }

    #[test]
    fn cone_row_matches_contact() {
        let gamma = SMatrix::<f64, N, M>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.0, 0.2, -0.3, 0.7, -0.4, 0.1,
            0.6,
        ]);
        let lam = SVector::<f64, N>::from_row_slice(&[0.5, -0.2, 0.8, 0.1, -0.6, 0.3]);
        let rows = Norm2.cone_constraints(&gamma);
        assert!(rows.linear.is_empty());
        assert_eq!(rows.soc.len(), 1);
        let (g, h) = &rows.soc[0];
        assert_relative_eq!(*h, 1.0, epsilon = 1e-15);
        // ||G lambda|| equals contact(Gamma^T lambda) equals the oracle value.
        assert_relative_eq!(
            (g * lam).norm(),
            Norm2.contact(gamma.transpose() * lam),
            epsilon = 1e-12
        );
        assert_relative_eq!((g * lam).norm(), 0.6428841264178172, epsilon = 1e-12);
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail (red)**

Run: `cargo test --lib cost::norm2`
Expected: FAIL — each test panics with `not implemented: Phase 2: g(y) = ||y||_2` (and the support/cone variants), because the stub bodies are still `unimplemented!`.

- [ ] **Step 3: Replace the stub implementation**

Replace everything **above** the `#[cfg(test)]` line you just added — the module doc comment, the `use` lines, the `Norm2` struct, and the `impl SublevelSet for Norm2` block — with:

```rust
//! L2 cost `||u||_2`: the unit-ball sublevel set (Table II).

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// L2 cost `||u||_2`. Contact `g(y) = ||y||_2`, support `s(y) = y / ||y||_2`,
/// and one SOC row per time: `||Gamma^T(t) lambda||_2 <= 1`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Norm2;

impl SublevelSet for Norm2 {
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        y.norm()
    }

    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        let n = y.norm();
        if n > 0.0 {
            y / n
        } else {
            SVector::<f64, M>::zeros()
        }
    }

    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        // g(Gamma^T lambda) <= 1  <=>  one SOC row with G = Gamma^T, h = 1.
        ConicRows {
            linear: Vec::new(),
            soc: vec![(gamma_t.transpose(), 1.0)],
        }
    }
}
```

- [ ] **Step 4: Run the tests to verify they pass (green)**

Run: `cargo test --lib cost::norm2`
Expected: PASS — `6 passed`.

- [ ] **Step 5: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
git add -A
git commit -m "feat(cost): implement Norm2 sublevel set (contact, support, SOC row)"
```

---

## Task 2: `FaceMax` — the tetrahedral face-max cost (`max(Vᶠᵃᶜᵉu)`)

Replace the stub in `src/cost/facemax.rs`. Adds a private `vertex_columns()` helper (the 4 `V_vertex` columns); `contact`/`support` range over `W = [0 | V_vertex]`, and `cone_constraints` emits the 4 linear rows.

**Files:**
- Modify: `src/cost/facemax.rs` (replace the stub `impl`; add `vertex_columns()` and a `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `super::SublevelSet`, `crate::types::{ConicRows, M, N}`, `nalgebra::{SMatrix, SVector}` (already imported in the stub).
- Produces:
  - `fn vertex_columns() -> [SVector<f64, M>; 4]` (module-private) — the 4 unit `V_vertex` columns.
  - `impl SublevelSet for FaceMax` with `contact` = `max(0, maxₖ y·vₖ)`, `support` = argmax column (lowest-index tie-break; origin/zero when all `y·vₖ ≤ 0`), `cone_constraints` = four linear rows `(Γ vₖ, 1.0)`.

- [ ] **Step 1: Append the failing tests to `src/cost/facemax.rs`**

Add this test module at the end of the stub file (the `vertex_columns()` helper it references will exist after Step 3; for now the stub bodies are `unimplemented!`, so these compile only once `vertex_columns` exists — so this step's "red" is the compile error in Step 2, then the panic disappears with Step 3). Add:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn vertices_are_unit_and_tetrahedral() {
        let cols = vertex_columns();
        for v in &cols {
            assert_relative_eq!(v.norm(), 1.0, epsilon = 1e-12);
        }
        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    assert_relative_eq!(cols[i].dot(&cols[j]), -1.0 / 3.0, epsilon = 1e-12);
                }
            }
        }
    }

    #[test]
    fn contact_known_directions() {
        let s23 = (2.0_f64 / 3.0).sqrt();
        let s13 = (1.0_f64 / 3.0).sqrt();
        assert_relative_eq!(
            FaceMax.contact(SVector::<f64, M>::new(1.0, 0.0, 0.0)),
            s23,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            FaceMax.contact(SVector::<f64, M>::new(0.0, 0.0, 1.0)),
            s13,
            epsilon = 1e-12
        );
    }

    #[test]
    fn contact_of_zero_is_zero() {
        assert_relative_eq!(
            FaceMax.contact(SVector::<f64, M>::zeros()),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn support_is_argmax_vertex() {
        let cols = vertex_columns();
        // y = (1,0,0): unique argmax = vertex column 0.
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::new(1.0, 0.0, 0.0)),
            cols[0],
            epsilon = 1e-12
        );
        // y = (0,0,1): tie between columns 2 and 3; lowest index (2) wins.
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::new(0.0, 0.0, 1.0)),
            cols[2],
            epsilon = 1e-12
        );
        // y = (-0.6,0.2,-0.9): argmax = vertex column 1.
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::new(-0.6, 0.2, -0.9)),
            cols[1],
            epsilon = 1e-12
        );
    }

    #[test]
    fn support_of_zero_is_zero() {
        assert_relative_eq!(
            FaceMax.support(SVector::<f64, M>::zeros()),
            SVector::<f64, M>::zeros(),
            epsilon = 1e-15
        );
    }

    #[test]
    fn contact_support_identity_eq23() {
        // eq. 23: lambda . s(lambda) = g(lambda).
        for y in [
            SVector::<f64, M>::new(0.3, 0.4, 0.5),
            SVector::<f64, M>::new(-0.6, 0.2, -0.9),
        ] {
            assert_relative_eq!(
                y.dot(&FaceMax.support(y)),
                FaceMax.contact(y),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn positive_homogeneity() {
        // g(a y) = a g(y) for a >= 0 (eq. 8 / Property 3).
        let y = SVector::<f64, M>::new(0.3, 0.4, 0.5);
        assert_relative_eq!(
            FaceMax.contact(y * 2.5),
            2.5 * FaceMax.contact(y),
            epsilon = 1e-12
        );
    }

    #[test]
    fn vertex_face_transcription_cross_check() {
        // V_face (eq. 48), 4x3. Used only to cross-check the transcription of
        // both matrices: f(v_k) = max_row(V_face v_k) is the same (= 1/9) for
        // every vertex column v_k of V_vertex.
        let s23 = (2.0_f64 / 3.0).sqrt();
        let s13 = (1.0_f64 / 3.0).sqrt();
        let v_face = SMatrix::<f64, 4, 3>::from_row_slice(&[
            -s23, 0.0, s13, s23, 0.0, s13, 0.0, -s23, -s13, 0.0, s23, -s13,
        ]) / 3.0;
        let cols = vertex_columns();
        for v in &cols {
            assert_relative_eq!((v_face * v).max(), 1.0 / 9.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn cone_rows_match_contact() {
        let gamma = SMatrix::<f64, N, M>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.0, 0.2, -0.3, 0.7, -0.4, 0.1,
            0.6,
        ]);
        let lam = SVector::<f64, N>::from_row_slice(&[0.5, -0.2, 0.8, 0.1, -0.6, 0.3]);
        let rows = FaceMax.cone_constraints(&gamma);
        assert!(rows.soc.is_empty());
        assert_eq!(rows.linear.len(), 4);
        for (_, b) in &rows.linear {
            assert_relative_eq!(*b, 1.0, epsilon = 1e-15);
        }
        let max_row = rows
            .linear
            .iter()
            .map(|(a, _)| a.dot(&lam))
            .fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(
            max_row,
            FaceMax.contact(gamma.transpose() * lam),
            epsilon = 1e-12
        );
        assert_relative_eq!(max_row, 0.37230594560185404, epsilon = 1e-12);
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail (red)**

Run: `cargo test --lib cost::facemax`
Expected: FAIL to compile — `cannot find function vertex_columns in this scope` (the helper does not exist yet). This is the intended red.

- [ ] **Step 3: Replace the stub implementation**

Replace everything **above** the `#[cfg(test)]` line — the module doc comment, the `use` lines, the `FaceMax` struct, and the `impl SublevelSet for FaceMax` block — with (note `vertex_columns()` is added between the struct and the impl):

```rust
//! Face-max cost `max(V_face u)` for the tetrahedral fixed-attitude occulter
//! (eq. 47-48). The algorithm consumes the cost only through its contact /
//! support / cone-constraint forms, all expressed via `W = [0 | V_vertex]`
//! (Table II); `V_face` itself is the cost definition and is not needed here.

use super::SublevelSet;
use crate::types::{ConicRows, M, N};
use nalgebra::{SMatrix, SVector};

/// Face-max cost. Contact `g(y) = max(0, max_k y . v_k)` over the four
/// `V_vertex` columns (the origin column of `W` supplies the `max(0, .)`);
/// support `s(y)` is the argmax column (origin -> zero when all `y . v_k <= 0`);
/// linear rows `(Gamma v_k)^T lambda <= 1` for each vertex column.
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceMax;

/// The four `V_vertex` columns (eq. 47): unit tetrahedral thruster directions.
fn vertex_columns() -> [SVector<f64, M>; 4] {
    let a = (2.0_f64 / 3.0).sqrt();
    let b = (1.0_f64 / 3.0).sqrt();
    [
        SVector::<f64, M>::new(a, 0.0, -b),
        SVector::<f64, M>::new(-a, 0.0, -b),
        SVector::<f64, M>::new(0.0, a, b),
        SVector::<f64, M>::new(0.0, -a, b),
    ]
}

impl SublevelSet for FaceMax {
    fn contact(&self, y: SVector<f64, M>) -> f64 {
        // g(y) = max over W = [0, V_vertex]; the origin column contributes 0.0,
        // so g(y) >= 0 always (and g(0) = 0).
        vertex_columns()
            .iter()
            .map(|v| y.dot(v))
            .fold(0.0, f64::max)
    }

    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M> {
        // argmax column of W; ties resolve to the lowest index; the origin
        // column (zero vector, value 0) wins when every y . v_k <= 0.
        let mut best = SVector::<f64, M>::zeros();
        let mut best_val = 0.0_f64;
        for v in vertex_columns() {
            let val = y.dot(&v);
            if val > best_val {
                best_val = val;
                best = v;
            }
        }
        best
    }

    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows {
        // g(Gamma^T lambda) <= 1  <=>  (Gamma v_k)^T lambda <= 1 for each vertex
        // column (the origin column gives the vacuous 0 <= 1 and is omitted).
        let linear = vertex_columns()
            .iter()
            .map(|v| (gamma_t * v, 1.0))
            .collect();
        ConicRows {
            linear,
            soc: Vec::new(),
        }
    }
}
```

- [ ] **Step 4: Run the tests to verify they pass (green)**

Run: `cargo test --lib cost::facemax`
Expected: PASS — `9 passed`.

- [ ] **Step 5: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
git add -A
git commit -m "feat(cost): implement FaceMax sublevel set (tetrahedral occulter, eq. 47-48)"
```

---

## Task 3: `Piecewise` — the eq. 49 time-varying selector

Replace the stub in `src/cost/piecewise.rs` and update the wiring test in `src/cost/mod.rs`. `Piecewise` holds a `Norm2`, a `FaceMax`, the orbit `period`, and the perigee-window `half_width`; `CostModel::at(t)` returns `&FaceMax` inside the perigee windows and `&Norm2` outside.

**Files:**
- Modify: `src/cost/piecewise.rs` (replace the stub struct + `impl`; add `#[cfg(test)] mod tests`)
- Modify: `src/cost/mod.rs` (update the `cost_types_wire_to_their_traits` test for `Piecewise::new`)

**Interfaces:**
- Consumes: `super::{CostModel, FaceMax, Norm2, SublevelSet}` (note the import list **grows** from the stub's `{CostModel, SublevelSet}`), plus `crate::types::M`, `nalgebra::SVector`, `approx` in tests.
- Produces:
  - `pub struct Piecewise { norm2: Norm2, facemax: FaceMax, period: f64, half_width: f64 }` (derives `Debug, Clone, Copy` — **not** `Default`).
  - `pub fn Piecewise::new(period: f64) -> Self` (sets `half_width = 3600.0`).
  - `pub fn Piecewise::in_perigee_window(&self, t: f64) -> bool` (eq. 49 set `T₁`).
  - `impl CostModel for Piecewise { fn at(&self, t: f64) -> &dyn SublevelSet }`.

- [ ] **Step 1: Fix the wiring test in `src/cost/mod.rs`**

`Piecewise` is about to gain fields, so the unit-struct form `&Piecewise` will no longer compile. Replace the test body:

```rust
    #[test]
    fn cost_types_wire_to_their_traits() {
        let _s: &dyn SublevelSet = &Norm2;
        let _f: &dyn SublevelSet = &FaceMax;
        // Piecewise now carries fields, so construct it via `new`.
        let pw = Piecewise::new(39_338.811_433_158_5);
        let _c: &dyn CostModel = &pw;
    }
```

- [ ] **Step 2: Append the failing tests to `src/cost/piecewise.rs`**

Add this test module at the end of the stub file:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::M;
    use approx::assert_relative_eq;
    use nalgebra::SVector;

    #[test]
    fn perigee_window_boundaries() {
        // period 40000 s: first perigee center at 20000, half_width 3600, so
        // T1 is the open interval (16400, 23600). Probe 1 s inside / outside
        // each edge. The exact 3600 s boundary is excluded by the strict `<`,
        // but it is a floating-point knife-edge, so it is not asserted directly.
        let pw = Piecewise::new(40000.0);
        // Center, and 3599 s either side -> inside T1.
        assert!(pw.in_perigee_window(20000.0));
        assert!(pw.in_perigee_window(16401.0));
        assert!(pw.in_perigee_window(23599.0));
        // 3601 s either side -> outside T1.
        assert!(!pw.in_perigee_window(16399.0));
        assert!(!pw.in_perigee_window(23601.0));
        // Apogees (t = 0 and t = period) -> outside T1.
        assert!(!pw.in_perigee_window(0.0));
        assert!(!pw.in_perigee_window(40000.0));
        // Second orbit's perigee at (1.5) * 40000 = 60000 -> inside T1.
        assert!(pw.in_perigee_window(60000.0));
    }

    #[test]
    fn at_selects_facemax_in_window_norm2_outside() {
        let pw = Piecewise::new(40000.0);
        let ex = SVector::<f64, M>::new(1.0, 0.0, 0.0);
        // Inside T1 -> FaceMax: g(ex) = sqrt(2/3).
        assert_relative_eq!(
            pw.at(20000.0).contact(ex),
            (2.0_f64 / 3.0).sqrt(),
            epsilon = 1e-12
        );
        // Outside T1 -> Norm2: g(ex) = 1.
        assert_relative_eq!(pw.at(0.0).contact(ex), 1.0, epsilon = 1e-12);
    }
}
```

- [ ] **Step 3: Run the tests to verify they fail (red)**

Run: `cargo test --lib cost::piecewise`
Expected: FAIL to compile — `no function or associated item named new found for struct Piecewise` / `no method named in_perigee_window` (the stub is a fieldless unit struct). The missing `new` surfaces from *both* the `cost/mod.rs` wiring edit (Step 1) and the appended `piecewise.rs` tests, so the whole `cost` module fails to compile — that is the intended red.

- [ ] **Step 4: Replace the stub implementation**

Replace everything **above** the `#[cfg(test)]` line — the module doc comment, the `use` line, the `Piecewise` struct, and the `impl CostModel for Piecewise` block — with:

```rust
//! Time-varying piecewise cost (eq. 49): FaceMax in 2-hr perigee windows (T1),
//! Norm2 elsewhere (T2).

use super::{CostModel, FaceMax, Norm2, SublevelSet};

/// Piecewise eq.-49 selector. `T1 = { t : |t - (k+0.5) period| < half_width }`
/// with `half_width = 1 hr` (eq. 49's 2-hr windows). The centers
/// `(k+0.5) period` land on perigee for the worked example because its chief
/// starts at apogee (`M0 = 180 deg`) at `t = 0`. The paper does not pin down
/// whether `period` is the Keplerian `2 pi / n` or the J2-perturbed period, so
/// the caller passes the period it wants (the worked example uses `2 pi / n`,
/// approx 10.93 hr, consistent with the paper's rounded 10.92 hr).
#[derive(Debug, Clone, Copy)]
pub struct Piecewise {
    norm2: Norm2,
    facemax: FaceMax,
    period: f64,
    half_width: f64,
}

impl Piecewise {
    /// Build the eq.-49 selector for an orbit period `period` [s]; the perigee
    /// window half-width is `1 hr = 3600 s`.
    pub fn new(period: f64) -> Self {
        Self {
            norm2: Norm2,
            facemax: FaceMax,
            period,
            half_width: 3600.0,
        }
    }

    /// `true` iff `t` lies within `half_width` of a perigee center
    /// `(k+0.5) period`, i.e. `t` is in the eq.-49 set `T1`.
    pub fn in_perigee_window(&self, t: f64) -> bool {
        let frac = (t / self.period).rem_euclid(1.0);
        (frac - 0.5).abs() * self.period < self.half_width
    }
}

impl CostModel for Piecewise {
    fn at(&self, t: f64) -> &dyn SublevelSet {
        if self.in_perigee_window(t) {
            &self.facemax
        } else {
            &self.norm2
        }
    }
}
```

- [ ] **Step 5: Run the tests to verify they pass (green)**

Run: `cargo test --lib cost`
Expected: PASS — the two `piecewise` tests, the updated `cost_types_wire_to_their_traits`, plus all Task 1/2 cost tests; `0 failed`.

- [ ] **Step 6: Run the whole suite to confirm nothing regressed**

Run: `cargo test`
Expected: PASS — all Phase 0/1 tests plus the Phase 2 cost tests; `0 failed`.

- [ ] **Step 7: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
git add -A
git commit -m "feat(cost): implement Piecewise eq.49 time-varying cost selector"
```

---

## Self-Review

**1. Spec coverage (design §5.2–5.3 + Phase 2 roadmap §6):**

| Spec requirement | Task |
|---|---|
| `Norm2` cost — `g(y)=‖y‖₂`, `s(y)=y/‖y‖₂` | Task 1 |
| `FaceMax` cost — `g(y)=maxₖ yᵀwₖ`, `s=argmax`, `W=[0,V_vertex]` | Task 2 |
| `V_vertex` / `V_face` (eq. 47–48) encoded | Task 2 (`vertex_columns`; `V_face` in the cross-check test) |
| `Piecewise` selector (eq. 49), `T₁/T₂` windows | Task 3 |
| Test: `λᵀs(λ)=g(λ)` (eq. 23) | Tasks 1 & 2 (`contact_support_identity_eq23`) |
| Test: positive homogeneity (eq. 8 / Property 3) | Tasks 1 & 2 (`positive_homogeneity`) |
| Test: known directions | Tasks 1 & 2 (`contact_*`/`support_*` known vectors) |
| Test: `T₁/T₂` window logic | Task 3 (`perigee_window_boundaries`, `at_selects_*`) |
| Table II "solver rows" (`cone_constraints`): 1 SOC / 4 linear | Tasks 1 & 2 (`cone_*_matches_contact`) — completes the trait, de-risks Phase 3 |
| Exit: cost tests pass | All tasks; verified green in a scratch worktree (18 cost tests) |

Every Phase 2 requirement maps to a task. `cone_constraints` is implemented now (the spec lists it under Table II "Solver rows" and the §4.1 trait comment calls it "what the solver layer needs"); it is pure cost math, independently tested against `contact`, and leaves no `unimplemented!` in the cost layer. The `V_face` matrix is realized only as a transcription cross-check (`f(vₖ)=1/9`) because the algorithm consumes the cost solely through `V_vertex`/`W`; this is the Rust↔spec self-consistency check, complementary to the char-by-char PDF read (done during planning).

**2. Placeholder scan:** No `TODO`/"fill in"/"handle edge cases". Every code step shows complete code; every reference number is concrete oracle output; every `Run:` has an expected result. Each cost's `unimplemented!`/`#[allow(unused_variables)]` stub is fully replaced.

**3. Type consistency:** `Norm2`, `FaceMax`, `Piecewise`, the free fn `vertex_columns() -> [SVector<f64,M>;4]`, `Piecewise::{new, in_perigee_window}`, and the trait methods `contact(&self, SVector<f64,M>) -> f64` / `support(&self, SVector<f64,M>) -> SVector<f64,M>` / `cone_constraints(&self, &SMatrix<f64,N,M>) -> ConicRows` / `at(&self, f64) -> &dyn SublevelSet` are named identically where defined and where consumed (the `Piecewise` tests and the `cost/mod.rs` wiring test). `ConicRows` is used exactly as defined in `types.rs` (`linear: Vec<(SVector<f64,N>, f64)>`, `soc: Vec<(SMatrix<f64,M,N>, f64)>`); Norm2 fills `soc` (one `(Γᵀ,1)`), FaceMax fills `linear` (four `(Γvₖ,1)`). Two stub changes are flagged explicitly rather than left to surprise the executor: `Piecewise` loses its `Default` derive and unit-struct form (Task 3 Step 1 updates the `cost/mod.rs` wiring test; Step 4 drops `Default`), and the `piecewise.rs` import list grows to include `FaceMax, Norm2`.
