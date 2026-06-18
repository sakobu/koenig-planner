# Phase 1 (J₂ mean-ROE dynamics) — Verification Report

- **Date:** 2026-06-17
- **Scope:** `src/dynamics/{constants,kepler,orbit,b_matrix,stm,j2_roe}.rs`
- **Verdict:** All `B(t)`, `Φ(t,t_f)`, eq.50 secular-rate, and constant terms verified correct. **Zero discrepancies** across five independent verification routes.

This report closes the Phase 1 exit gate from the design spec: *"every STM/`B` term verified character-by-character against `docs/Planner.pdf`."*

## Verification routes

| # | Route | Result |
|---|---|---|
| 1 | Entrywise vs a Python oracle transcribed term-for-term from spec §5.4 | Rust ↔ Python agree ≤ 1e-9 (`B`, `Φ`, `Γ=Φ·B`) |
| 2 | Character-by-character read of every term vs `docs/Planner.pdf` p.12–13 (by hand) | All match |
| 3 | 10-agent adversarial re-read vs `Planner.pdf` (redundant) + primary-source triangulation | All `all_match`, 0 mismatch, 0 uncertain |
| 4 | Independent finite-difference of `B(t)` via Cartesian r,v (uses none of `B`'s formulas) | ≤ 1e-9 Frobenius (`tests/fd_b_matrix.rs`) |
| 5 | Convention/limit tests | `Φ→I` as `Δt→0`; `Φ₂₄`≠0; `√(a/μ)` scaling; Kepler round-trip + known-`ν` pairs; secular/`n` anchors |

## Primary-source triangulation (route 3)

- **`B(t)` ⟷ Chernick & D'Amico 2018, Eq. (38), p.23.** Prefactor `1/(na) = √(a/μ)` confirmed; all 11 nonzero terms and the zero structure match.
- **`Φ(t,t_f)` ⟷ Koenig–Guffanti–D'Amico 2017 (ref [27]), Appendix A6.** Rows 1, 3, 4, 5, 6 identical to A6. Row 2 (δλ) is the **documented Koenig-2020 modification**: ref [27] A6 prints `(7/2)κ·E·P` with `E=1+η`; Koenig 2020 swaps in the η-weighted δλ giving `Φ₂₁=(−1.5nΔt−7κηP)Δt`, `Φ₂₃=7κe_{x1}PΔt/η`, `Φ₂₄=7κe_{y1}PΔt/η`, `Φ₂₅=−7κηSΔt` — matched verbatim by the code.

## Subtleties independently confirmed (the easy-to-get-wrong spots)

- **`Ṁ` (eq.50) uses `η³`, not `η⁴`** in the J₂ term — flagged and confirmed.
- **`Φ` mixed initial/final eccentricity subscripts:** `Φ₃₃` uses `e_{x1}·e_{y2}`, `Φ₄₄` uses `e_{y1}·e_{x2}` (asymmetric); `Φ₂₃/Φ₂₄/Φ₆₃/Φ₆₄` use initial (`x1/y1`); `Φ₃₁/Φ₄₁/Φ₃₅/Φ₄₅` use final (`x2/y2`). Code matches the paper's exact assignment.
- **`B₃₃` positive vs `B₄₃` negative** sign asymmetry preserved.
- **Mean motion** is `√(μ/a³)` (p.13 `Φ₂₁` prints it explicitly), resolving any `a²/a³` rendering ambiguity in eq.50.
- **`Φ₂₄` is nonzero** (`7κe_{y1}PΔt/η`): δλ couples to δe_y; matched in code and confirmed in the printed term list.

## Bug caught and fixed during verification

At `M=π` (apoapsis), `wrap_to_pi(π) = −π`, so the Kepler solver returns `ν = −π` (physically `±π`, same point). The Kepler known-pairs test was corrected to compare `|ν|`. Caught by the pre-execution scratch run; no impact on `B`/`Φ`/`Γ`.

## Test inventory (27 green)

- `dynamics::kepler` (3): wrap, known `M→ν` pairs, Kepler residual round-trip.
- `dynamics::orbit` (3): `n`/`η` anchors, secular-rate anchors, propagation linearity.
- `dynamics::b_matrix` (3): zero structure, entrywise oracle, `√(a/μ)` scaling law.
- `dynamics::stm` (3): `Φ→I`, `Φ₂₄`≠0, entrywise oracle.
- `dynamics::j2_roe` (3): trait-object, `Γ(t_f)=B(t_f)`, entrywise `Γ` oracle.
- `tests/fd_b_matrix.rs` (1): independent finite-difference `B`.
- Plus 11 Phase 0 type/API tests.
