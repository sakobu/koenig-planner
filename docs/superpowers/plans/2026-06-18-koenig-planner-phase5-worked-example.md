# Phase 5 — Worked-Example Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the paper's published worked example (Table III inputs → Table IV's 3-maneuver, ≈82.4 mm/s solution) end-to-end through the existing `solve(...)` pipeline, as a CI-enforced integration test plus a human-facing `examples/mdot.rs` report and Fig. 7 contact-curve data — and close the two Phase-4-deferred refinement-test items.

**Architecture:** Phase 5 writes **no new library code** for the primary goal — `solve` already exists and was whole-branch-reviewed in Phase 4. The work is (1) *encoding* the Table III inputs with the correct `a_c` scaling, (2) *running* the full pipeline against the real ill-conditioned `J2Roe` `Γ` for the first time, and (3) *asserting* the §7 targets within solver-tolerance bands. Because this is the first time the whole pipeline meets published numbers, **every numeric task is characterization-first**: run, print, compare to §7, and only then lock bands — a mismatch is a debugging signal (systematic-debugging skill), never a reason to silently widen a band. Two small `src/algorithm/refine.rs` changes close the Phase-4 deferrals.

**Tech Stack:** Rust 2021 (rust-version 1.92), `nalgebra 0.35` (`SVector`/`SMatrix`), `clarabel 0.11` (conic solver, used transitively via `solve`), `csv`/`plotters` (behind the `validation` feature, Fig. 7 output only), `approx 0.5` (dev, test assertions). CI = GitHub Actions: `fmt` + `clippy --all-features -D warnings` + `build --all-features` + `test --all-features`.

## Global Constraints

- **Work on a feature branch, not `main`.** Branch name: `phase5-worked-example`. Prior phases each merged via a PR (Phase 1–4 = PRs #1/#3/#5). Open a PR to `main` at the end; the CI gate (`fmt` + `clippy --all-features -D warnings` + `build --all-features` + `test --all-features`) must be green before merge. Close out with the `superpowers:finishing-a-development-branch` skill.
- **Source of truth for numbers:** `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md` §7 (worked-example targets) and §5.5 (units & `a_c` scaling). The PDF is `docs/Planner.pdf`.
- **Bands, not bit-equality.** Spec §7 + Risk 3/8: "Assertions use sensible numerical bands … not bit-for-bit equality — exact figures depend on solver tolerances." Every numeric assertion is a band tied to the paper's `ε = 0.01` tolerances. The starting bands in this plan are §7-derived; a task may *widen* a band only with a one-line comment justifying it as solver tolerance (cite the observed value), and may **never** widen a band to mask a structural discrepancy (wrong maneuver count, wrong times, residual far above target) — that is a systematic-debugging case.
- **`a_c` scaling (spec §5.5 — the highest-risk decision):** `a_c = 25_000e3 m` (25 000 km). The native pseudostate the solver reconstructs is **dimensionless** ROE: `w_nd = w_metres / a_c`. Table III lists `w` in metres only because it is `a_c · w_nd` for display. **Feed `w_nd` (dimensionless) to `solve`.** Maneuver Δv then comes out in **m/s** directly (report ×1000 for mm/s). Do **not** bake `a_c` into `B`/`Φ`/`Γ`. Consequence for the dual: §7's `λ_opt ≈ 1e-6·[…]` pairs with the *metre-scaled* `w`; solving with `w_nd` yields `λ` scaled by `a_c` (≈ `25·[…]`). Only `λ`'s **direction** and the ratio `λᵀw/g(λ)` are scale-invariant — assert on direction.
- **Ill-conditioning (the central Phase-5 risk):** `Γ(t)`'s δλ row (row index 1) is ~`1e3` while every other row is ~`1e-4` (condition number ~`1e7`; verified: `J2Roe::gamma(16_050)` row 1 = `[1068.6, -1152.8, 2.1e-6]`). The optimal maneuvers' large δλ contributions **nearly cancel** to reach the tiny target δλ. clarabel's default Ruiz equilibration is expected to handle the row-scale disparity, but residual/`total_dv` landing slightly off the paper is a *tolerance* question, not necessarily a bug — characterize before judging.
- **No `cargo test` of examples.** CI `test` runs the `tests/` integration files (authoritative gate); CI `build`/`clippy` compile-and-lint `examples/mdot.rs` (incl. its `validation`-gated `csv` code) but do **not** run it. So the §7 *gate* lives in `tests/worked_example.rs`; the example's own `assert!`s are belt-and-suspenders, fired only by manual `cargo run --example mdot`.
- **`approx` is a dev-dependency** → usable in `tests/`, **not** in `examples/`. The example uses plain float comparisons.

---

## File Structure

| File | Disposition | Responsibility |
|---|---|---|
| `tests/worked_example.rs` | **Create** | CI-enforced §7 assertions: the primary Table III/IV case (Task 1–2) and the secondary Hunter cross-check (Task 4). |
| `examples/mdot.rs` | **Replace** the 9-line Phase-0 stub | Human-facing Table-IV report + Fig. 7 contact-function curve (`(t, g)` CSV behind `--features validation`). Self-checks §7 via `assert!`. (Task 3) |
| `src/algorithm/refine.rs` | **Modify** (`RefineOutcome` + its `#[cfg(test)]` module) | Add `active_set_trace` observability field (Task 5); add the deferred real-`J2Roe` ≥3-iteration drop-then-readd test (Task 5) and the `achieved > target` assertion (Task 6). |
| `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md` | **Modify** | Mark Phase 5 done; record actual numbers + scaling/tolerance findings (Task 7). |
| `.claude/.../memory/` | **Modify** | Update `spec-validation-status.md` + `MEMORY.md` index (Task 7). |

**Reused public API (verified on disk — do not re-derive):**

```rust
// crate root re-exports (src/lib.rs)
koenig_planner::{solve, Solution, SolveParams, TimeGrid, Maneuver, Pseudostate, Dual,
                 PlannerError, CostModel, SublevelSet, Dynamics, N, M};
koenig_planner::dynamics::{J2Roe, AbsoluteOrbit};   // not re-exported at crate root
koenig_planner::cost::Piecewise;                    // not re-exported at crate root

// solve (src/algorithm/mod.rs:56) — feed w as the DIMENSIONLESS w_nd
pub fn solve<D: Dynamics, C: CostModel>(
    dynamics: &D, cost: &C, w: Pseudostate, grid: TimeGrid, params: &SolveParams,
) -> Result<Solution, PlannerError>;

// Solution (src/types.rs:98)
pub struct Solution {
    pub maneuvers: Vec<Maneuver>, // sorted by t ascending (t_opt is sorted before extract)
    pub total_dv: f64,            // Σ‖Δvⱼ‖ [m/s]
    pub iterations: usize,        // Algorithm 2 solve count
    pub residual: f64,            // ‖w_err‖/‖w‖ (dimensionless ratio — scale-invariant)
    pub lambda: Dual,             // SVector<f64,6>; direction is the scale-invariant part
}
pub struct Maneuver { pub t: f64 /* s */, pub dv: SVector<f64, 3> /* m/s, RTN */ }

// AbsoluteOrbit::new(a[m], e, i[rad], raan[rad], argp[rad], mean_anom[rad])  (orbit.rs:38)
// AbsoluteOrbit::mean_motion(&self) -> f64   // n = sqrt(mu/a^3) [rad/s]      (orbit.rs:50)  [pub]
// J2Roe::new(chief_ti: AbsoluteOrbit, t_i: f64, t_f: f64) -> J2Roe           (j2_roe.rs:23)
// Piecewise::new(period: f64 /* s */) -> Piecewise                           (piecewise.rs:24)
// TimeGrid::uniform(t_i, t_f, dt) ; .len() ; .time(idx) ; .times()           (types.rs:43)
// SolveParams::default() == Table III: n_coarse=20, n_init=6, eps_*=0.01, q=I (types.rs:85)
```

---

## Task 1: Worked-example integration test — characterization (no §7 bands yet)

Stand up `tests/worked_example.rs`, encode Table III, run the **real** pipeline once, and **print** the full solution. This de-risks the `a_c` scaling and ill-conditioning *before* any band is committed. The committed deliverable asserts only that the solve succeeds and is finite; Task 2 calibrates the §7 bands against what this prints.

**Files:**
- Create: `tests/worked_example.rs`

**Interfaces:**
- Consumes: `solve`, `Solution`, `SolveParams`, `TimeGrid`, `Pseudostate` (crate root); `J2Roe`, `AbsoluteOrbit` (`dynamics::`); `Piecewise` (`cost::`).
- Produces: a module-private `fn worked_example_inputs() -> (J2Roe, Piecewise, Pseudostate, TimeGrid, SolveParams)` and `const A_C: f64` reused by Task 2's assertions. (The `examples/mdot.rs` copy in Task 3 duplicates this small block — Rust makes sharing a helper between a `tests/` file and an `examples/` binary awkward, and the block is ~12 lines; duplication is the YAGNI choice here.)

- [ ] **Step 1: Create the branch**

```bash
git checkout -b phase5-worked-example
```

- [ ] **Step 2: Write the characterization test (prints everything, asserts only finiteness)**

Create `tests/worked_example.rs`:

```rust
//! Phase 5 worked-example validation: Table III inputs -> Table IV outputs
//! (paper §7). Characterization-first: Task 1 prints the real solution; Task 2
//! locks the §7 bands. Run with `--nocapture` to see the printout.

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_planner::{solve, Pseudostate, SolveParams, TimeGrid};
use std::f64::consts::TAU;

/// Chief semimajor axis a_c [m] — the I/O scaling factor (spec §5.5).
const A_C: f64 = 25_000e3;

/// Table III target pseudostate, displayed in metres (= a_c * w_nd).
const W_METRES: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

/// Table III inputs, ready to hand to `solve`. `w` is the DIMENSIONLESS w_nd.
fn worked_example_inputs() -> (J2Roe, Piecewise, Pseudostate, TimeGrid, SolveParams) {
    // Mean absolute orbit (Table III): a=25000 km, e=0.7, i=40°, Ω=358°, ω=0°, M=180°.
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let t_i = 0.0;
    let t_f = 117_990.0; // 3 orbits
    let dynamics = J2Roe::new(chief, t_i, t_f);
    // eq.49 perigee windows keyed to the Keplerian period 2π/n (≈ 39338.8 s).
    let cost = Piecewise::new(TAU / chief.mean_motion());
    // Feed the dimensionless w_nd = w_metres / a_c (spec §5.5).
    let w = Pseudostate::from_row_slice(&W_METRES) / A_C;
    let grid = TimeGrid::uniform(t_i, t_f, 30.0); // 3934 candidate times
    let params = SolveParams::default(); // Table III: n_coarse=20, n_init=6, eps=0.01, Q=I
    (dynamics, cost, w, grid, params)
}

#[test]
fn worked_example_characterization() {
    let (dynamics, cost, w, grid, params) = worked_example_inputs();
    assert_eq!(grid.len(), 3934, "Table III grid should have 3934 times");

    let sol = solve(&dynamics, &cost, w, grid, &params).expect("worked example should solve");

    eprintln!("=== Phase 5 worked-example characterization ===");
    eprintln!("iterations = {}", sol.iterations);
    eprintln!("total_dv   = {:.6} m/s  ({:.4} mm/s)", sol.total_dv, sol.total_dv * 1e3);
    eprintln!("residual   = {:.3e}", sol.residual);
    eprintln!("maneuvers  = {}", sol.maneuvers.len());
    for (j, m) in sol.maneuvers.iter().enumerate() {
        eprintln!(
            "  [{j}] t = {:>9.1} s   u = [{:>8.3}, {:>8.3}, {:>8.3}] mm/s   |Δv| = {:.4} mm/s",
            m.t, m.dv[0] * 1e3, m.dv[1] * 1e3, m.dv[2] * 1e3, m.dv.norm() * 1e3,
        );
    }
    eprintln!("lambda          = {:?}", sol.lambda.as_slice());
    eprintln!("lambda / a_c    = {:?}", (sol.lambda / A_C).as_slice());

    // Characterization only: assert success + finiteness, NOT §7 bands (Task 2).
    assert!(!sol.maneuvers.is_empty());
    assert!(sol.total_dv.is_finite() && sol.total_dv > 0.0);
    assert!(sol.residual.is_finite());
    assert!(sol.lambda.iter().all(|x| x.is_finite()));
}
```

- [ ] **Step 3: Run it and FAIL fast if it does not even solve**

Run: `cargo test --test worked_example worked_example_characterization -- --nocapture`
Expected: PASS, and a printed block. **Record the printed numbers** (paste them into the eventual commit message / Task 7 notes).

- [ ] **Step 4: Compare the printout to spec §7 (decision point — do not skip)**

§7 expects: 3 maneuvers at `t = [16050, 23280, 107100] s`; `u_R=[9.68,0.00,16.51]`, `u_T=[−23.02,−0.40,15.68]`, `u_N=[−25.56,−0.04,40.26]` mm/s; `total_dv ≈ 82.4 mm/s`; `~3 iterations`; `residual < 0.01%`; `λ_opt ∝ [34.97,3.42,30.68,17.84,−9.34,146.79]`.

- **If the printout matches §7 within ~1–2%** (3 maneuvers, times within a grid cell, `total_dv ≈ 0.082`): scaling and conditioning are confirmed → proceed to Task 2.
- **If `total_dv` (or the maneuver Δv) is off by a clean factor of `a_c` (= 25e6) or its reciprocal:** the scaling is inverted — re-read spec §5.5, confirm `w = w_metres / a_c`, and that the maneuver `dv` is reported as `dv * 1e3` for mm/s (not `dv / 1e3`). Fix and re-run Step 3.
- **If the structure is wrong** (≠3 maneuvers, times nowhere near §7, residual ≫ `1e-3`, or `solve` returns `NotConverged`): **STOP and invoke `superpowers:systematic-debugging`.** This is a real discrepancy in the pipeline against published numbers — likely candidates: the `Piecewise` period, a dynamics/scaling convention, or solver conditioning. Do **not** proceed to Task 2 or weaken anything until the root cause is understood.

- [ ] **Step 5: Commit the characterization**

```bash
git add tests/worked_example.rs
git commit -m "test(phase5): worked-example characterization (Table III -> solve, prints Table IV)"
```

---

## Task 2: Lock the §7 worked-example assertions

Calibrate the §7 bands against Task 1's printout and add the authoritative, CI-enforced assertions to the same test file.

**Files:**
- Modify: `tests/worked_example.rs`

**Interfaces:**
- Consumes: `worked_example_inputs()` + `A_C` from Task 1.
- Produces: the `worked_example_matches_table_iv` test (the Phase-5 exit gate).

- [ ] **Step 1: Add the §7-band assertion test**

Append to `tests/worked_example.rs` (add `use nalgebra::SVector;` to the imports):

```rust
#[test]
fn worked_example_matches_table_iv() {
    let (dynamics, cost, w, grid, params) = worked_example_inputs();
    let sol = solve(&dynamics, &cost, w, grid, &params).expect("worked example should solve");

    // §7: exactly 3 maneuvers, ascending in t. (solve sorts t_opt before extract.)
    assert_eq!(sol.maneuvers.len(), 3, "expected 3 maneuvers, got {}", sol.maneuvers.len());

    // §7 Table IV, mm/s, in t-ascending order: (t_s, [u_R, u_T, u_N]).
    let expected: [(f64, [f64; 3]); 3] = [
        (16_050.0, [9.68, -23.02, -25.56]),
        (23_280.0, [0.00, -0.40, -0.04]),
        (107_100.0, [16.51, 15.68, 40.26]),
    ];
    for (m, (t_exp, dv_mm_exp)) in sol.maneuvers.iter().zip(expected.iter()) {
        // Times are grid points; allow ±1 grid cell (30 s) of solver slack.
        assert!((m.t - t_exp).abs() <= 30.0, "maneuver time {} s vs §7 {} s", m.t, t_exp);
        for c in 0..3 {
            let got_mm = m.dv[c] * 1e3;
            // ±0.8 mm/s absolute band (~1–2% of the larger components; loosely
            // bounds the near-zero second maneuver). Widen ONLY with a documented
            // solver-tolerance reason (cite the observed value); never to hide a bug.
            assert!(
                (got_mm - dv_mm_exp[c]).abs() <= 0.8,
                "maneuver at {:.0}s component {c}: {got_mm:.3} mm/s vs §7 {:.2} mm/s",
                m.t, dv_mm_exp[c],
            );
        }
    }

    // §7: total Δv ≈ 82.4 mm/s; dual lower bound 82.0; "≤ 1% above" ⇒ ≤ 82.82.
    let total_mm = sol.total_dv * 1e3;
    assert!(
        (81.5..=82.8).contains(&total_mm),
        "total Δv = {total_mm:.4} mm/s, expected ≈ 82.4 (band 81.5–82.8)",
    );

    // §7: "~3 iterations" of Algorithm 2 (band per Risk 3).
    assert!((2..=6).contains(&sol.iterations), "iterations = {} (expected ~3)", sol.iterations);

    // §7: residual < 0.01% = 1e-4. Band 5e-4 absorbs clarabel's budget-binding
    // tolerance (Phase-4 note b); >5e-4 is a triage signal, not a band to widen.
    assert!(sol.residual < 5e-4, "residual = {:.3e}, §7 target < 1e-4", sol.residual);

    // §7: λ_opt ∝ [34.97, 3.42, 30.68, 17.84, −9.34, 146.79]. Solving with the
    // dimensionless w_nd scales λ by a_c, so compare DIRECTION (scale-invariant).
    let lambda_ref = SVector::<f64, 6>::from_row_slice(&[
        34.97, 3.42, 30.68, 17.84, -9.34, 146.79,
    ]);
    let cos = sol.lambda.dot(&lambda_ref) / (sol.lambda.norm() * lambda_ref.norm());
    assert!(cos > 0.999, "λ direction cosine vs §7 = {cos:.6} (expected > 0.999)");
}
```

- [ ] **Step 2: Run the gate test**

Run: `cargo test --test worked_example worked_example_matches_table_iv -- --nocapture`
Expected: PASS. If any band fails, apply the Step-4 decision rules from Task 1 (tolerance-widen with a documented reason, or systematic-debugging for a structural miss).

- [ ] **Step 3: Run the whole file + clippy on it**

Run: `cargo test --test worked_example`
Expected: PASS (both tests).
Run: `cargo clippy --all-features --tests -- -D warnings`
Expected: clean (no warnings).

- [ ] **Step 4: Commit**

```bash
git add tests/worked_example.rs
git commit -m "test(phase5): assert Table IV worked-example targets within §7 bands"
```

---

## Task 3: `examples/mdot.rs` — Table IV report + Fig. 7 contact curve

Replace the Phase-0 stub with the runnable worked example: a formatted Table-IV report, the Fig. 7 contact-function curve `g(t)` over the grid, CSV output behind `--features validation`, and self-checking `assert!`s.

**Files:**
- Replace: `examples/mdot.rs` (currently a 9-line stub)

**Interfaces:**
- Consumes: `solve`, `Solution`, `SolveParams`, `TimeGrid`, `Pseudostate` (root); `J2Roe`, `AbsoluteOrbit`, `Dynamics` (`dynamics::`); `Piecewise`, `CostModel`, `SublevelSet` (`cost::`/root). `Dynamics::gamma`, `CostModel::at`, `SublevelSet::contact` are needed for the Fig. 7 curve, so their traits must be in scope.
- Produces: a binary that prints the report and (with `--features validation`) writes `target/fig7_contact.csv`. No code consumed by later tasks.

- [ ] **Step 1: Replace the stub**

Overwrite `examples/mdot.rs`:

```rust
//! Worked example (paper §7): Table III inputs -> Table IV maneuvers, plus the
//! Fig. 7 contact-function curve. Run:
//!
//!   cargo run --example mdot                      # report only
//!   cargo run --example mdot --features validation # also writes target/fig7_contact.csv
//!
//! The CI-enforced §7 assertions live in tests/worked_example.rs; the asserts
//! here make a manual run self-checking.

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_planner::{solve, CostModel, Dynamics, Pseudostate, SolveParams, SublevelSet, TimeGrid};
use std::f64::consts::TAU;

/// Chief semimajor axis a_c [m] — the I/O scaling factor (spec §5.5).
const A_C: f64 = 25_000e3;
/// Table III target pseudostate in metres (= a_c * w_nd).
const W_METRES: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

fn main() {
    // --- Table III inputs (duplicated from tests/worked_example.rs; see plan note). ---
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let (t_i, t_f) = (0.0, 117_990.0);
    let dynamics = J2Roe::new(chief, t_i, t_f);
    let cost = Piecewise::new(TAU / chief.mean_motion());
    let w = Pseudostate::from_row_slice(&W_METRES) / A_C;
    let grid = TimeGrid::uniform(t_i, t_f, 30.0);
    let params = SolveParams::default();

    let sol = solve(&dynamics, &cost, w, grid, &params).expect("worked example should solve");

    // --- Table IV report. ---
    println!("Koenig planner — worked example (Table III -> Table IV)");
    println!("  candidate times : {}", grid.len());
    println!("  iterations      : {}", sol.iterations);
    println!("  residual        : {:.3e}", sol.residual);
    println!("  total Δv        : {:.4} mm/s", sol.total_dv * 1e3);
    println!("  maneuvers ({}):", sol.maneuvers.len());
    println!("    {:>10}  {:>9}  {:>9}  {:>9}  {:>9}", "t [s]", "u_R", "u_T", "u_N", "|Δv|");
    println!("    {:>10}  {:>9}  {:>9}  {:>9}  {:>9}", "", "[mm/s]", "[mm/s]", "[mm/s]", "[mm/s]");
    for m in &sol.maneuvers {
        println!(
            "    {:>10.0}  {:>9.2}  {:>9.2}  {:>9.2}  {:>9.3}",
            m.t, m.dv[0] * 1e3, m.dv[1] * 1e3, m.dv[2] * 1e3, m.dv.norm() * 1e3,
        );
    }
    println!("  λ_opt           : {:?}", sol.lambda.as_slice());
    println!("  λ_opt / a_c     : {:?}", (sol.lambda / A_C).as_slice());

    // --- Fig. 7: contact function g(t) = g_{U(1,t)}(Γᵀ(t) λ_opt) over the grid. ---
    let curve: Vec<(f64, f64)> = grid
        .times()
        .map(|t| {
            let y = dynamics.gamma(t).transpose() * sol.lambda;
            (t, cost.at(t).contact(y))
        })
        .collect();
    let max_g = curve.iter().map(|&(_, g)| g).fold(f64::NEG_INFINITY, f64::max);
    println!("  Fig. 7 max_t g  : {:.6}  (should be ≤ 1 + eps_cost = 1.01)", max_g);

    #[cfg(feature = "validation")]
    {
        let path = "target/fig7_contact.csv";
        let mut wtr = csv::Writer::from_path(path).expect("open fig7 csv");
        wtr.write_record(["t_s", "g"]).expect("write header");
        for (t, g) in &curve {
            wtr.write_record(&[t.to_string(), g.to_string()]).expect("write row");
        }
        wtr.flush().expect("flush fig7 csv");
        println!("  Fig. 7 curve    : written to {path} ({} rows)", curve.len());
    }
    #[cfg(not(feature = "validation"))]
    println!("  Fig. 7 curve    : (build with --features validation to write target/fig7_contact.csv)");

    // --- Self-check the headline §7 targets (belt-and-suspenders; gate is the test). ---
    assert_eq!(sol.maneuvers.len(), 3, "expected 3 maneuvers");
    assert!((81.5..=82.8).contains(&(sol.total_dv * 1e3)), "total Δv out of §7 band");
    assert!(sol.residual < 5e-4, "residual above §7 band");
    assert!(max_g <= 1.0 + params.eps_cost + 1e-9, "Fig. 7 max_t g exceeds 1 + eps_cost");
}
```

- [ ] **Step 2: Run without the feature**

Run: `cargo run --example mdot`
Expected: the report prints, the final line notes the `--features validation` hint, and the process exits 0 (all `assert!`s pass). The printed maneuvers/total should match the §7 numbers seen in Task 1.

- [ ] **Step 3: Run with the validation feature**

Run: `cargo run --example mdot --features validation`
Expected: same report plus `Fig. 7 curve : written to target/fig7_contact.csv (3934 rows)`.
Verify: `head -3 target/fig7_contact.csv` shows a `t_s,g` header and `(t, g)` rows; spot-check that `g` peaks near `1.0` at the three maneuver times (e.g. `grep "^16050," target/fig7_contact.csv`). (`target/` is already git-ignored — the CSV is not committed.)

- [ ] **Step 4: Lint the example under all features**

Run: `cargo clippy --all-features --examples -- -D warnings`
Expected: clean. Then `cargo build --all-features` (mirrors CI; compiles the `validation`-gated `csv` path).
Expected: builds clean.

- [ ] **Step 5: Commit**

```bash
git add examples/mdot.rs
git commit -m "feat(phase5): examples/mdot worked-example report + Fig. 7 contact curve"
```

---

## Task 4: Hunter & D'Amico second worked example (secondary cross-check)

Add the independent integration case from spec §7 ("Second worked example"). It uses the *same* J₂ ROE dynamics with a different orbit/target, so it is additional coverage; the Table III/IV case stays primary. The chief is given in `(e_x, e_y, u₀)` form and must be converted to `[e, ω, M]` — flagged as the one ambiguous step, with the dual lower bound as the decisive check.

**Files:**
- Modify: `tests/worked_example.rs`

**Interfaces:**
- Consumes: same imports as Task 1/2 (already in the file).
- Produces: the `hunter_cross_check` test.

- [ ] **Step 1: Add the Hunter test (characterization-aware)**

Append to `tests/worked_example.rs`:

```rust
#[test]
fn hunter_cross_check() {
    // Hunter & D'Amico 2025 "Sequential Formulation Validation" (spec §7).
    // Chief given as e_x, e_y, u₀ ⇒ convert: e = ‖(e_x,e_y)‖, ω = atan2(e_y,e_x),
    // M = u₀ − ω (u₀ = mean argument of latitude = ω + M). If the output misses,
    // the alternative reading is u₀ = TRUE arg of latitude (ω + ν) — try that next.
    let (e_x, e_y) = (-0.658, -0.239);
    let e = (e_x * e_x + e_y * e_y).sqrt(); // ≈ 0.7001
    let argp = e_y.atan2(e_x); // ≈ 200° (−2.793 rad)
    let u0 = 65.0_f64.to_radians();
    let mean_anom = u0 - argp; // M = u₀ − ω
    let chief = AbsoluteOrbit::new(
        A_C, // 25 000 km
        e,
        51.0_f64.to_radians(),
        30.0_f64.to_radians(),
        argp,
        mean_anom,
    );
    let (t_i, t_f) = (0.0, 39_000.0); // ~1 orbit
    let dynamics = J2Roe::new(chief, t_i, t_f);
    let cost = Piecewise::new(TAU / chief.mean_motion());
    let w = Pseudostate::from_row_slice(&[0.66, -1.52, -0.38, -1.44, 0.29, -0.91]) / A_C;
    let grid = TimeGrid::uniform(t_i, t_f, 10.0); // 3901 candidate times
    assert_eq!(grid.len(), 3901);
    let params = SolveParams::default();

    let sol = solve(&dynamics, &cost, w, grid, &params).expect("Hunter case should solve");
    eprintln!(
        "[hunter] iters={} total_dv={:.4e} m/s residual={:.3e} maneuvers={}",
        sol.iterations, sol.total_dv, sol.residual, sol.maneuvers.len(),
    );

    // §7: 3 maneuvers, total Δv ≈ 23.03e-5 m/s, dual lower bound 22.94e-5, ~4 iters.
    assert_eq!(sol.maneuvers.len(), 3, "expected 3 maneuvers, got {}", sol.maneuvers.len());
    assert!(
        (2.28e-4..=2.33e-4).contains(&sol.total_dv),
        "total Δv = {:.4e} m/s, expected ≈ 2.303e-4 (≥ dual bound 2.294e-4)",
        sol.total_dv,
    );
    assert!((2..=7).contains(&sol.iterations), "iterations = {} (expected ~4)", sol.iterations);
    assert!(sol.residual < 5e-4, "residual = {:.3e}, §7 target < 1e-4", sol.residual);
}
```

- [ ] **Step 2: Run it (characterize first if it fails)**

Run: `cargo test --test worked_example hunter_cross_check -- --nocapture`
Expected: PASS. If the maneuver count or `total_dv` misses while the primary case (Tasks 1–2) passes, the angle conversion is the suspect: switch `u₀` to the true-anomaly reading (compute `ν` from `M` is not needed here — instead set `mean_anom` so that `ω + M` matches; pragmatically, try `let mean_anom = u0 - argp;` vs treating `u0` as true latitude and back-solving). Use the `eprintln!` line + the dual lower bound (`total_dv ≥ 2.294e-4`) as the oracle. Mark with a comment whichever reading reproduces the paper.

- [ ] **Step 3: Commit**

```bash
git add tests/worked_example.rs
git commit -m "test(phase5): Hunter & D'Amico second worked example (secondary cross-check)"
```

---

## Task 5: Deferred Phase-4 item (a) — real-`J2Roe` ≥3-iteration drop-then-readd refine test

Phase 4's whole-branch review deferred to Phase 5 "a refinement test on the real ill-conditioned `J2Roe` `Γ` that observably runs ≥3 iterations with a drop-then-readd." `RefineOutcome` currently exposes only `max_g_trace` and the final `t_opt`, so the per-iteration active-set sizes are not observable. Add a small `active_set_trace` field, then write the test.

**Files:**
- Modify: `src/algorithm/refine.rs` (`RefineOutcome` struct + the `refine` loop + the `#[cfg(test)]` module)

**Interfaces:**
- Consumes (test-internal): `refine`, `RefineOutcome`, `J2Roe`, `AbsoluteOrbit`, `TimeGrid`, `Piecewise`, `SolveParams`. The test lives **inside** `src/algorithm/refine.rs` because `refine`/`RefineOutcome` are `pub(super)` (not reachable from `tests/`).
- Produces: `RefineOutcome.active_set_trace: Vec<usize>` (size of `T^est` solved at each iteration), `#[allow(dead_code)]` (read only by tests, mirroring `max_g_trace`).

- [ ] **Step 1: Add the `active_set_trace` field to `RefineOutcome`**

In `src/algorithm/refine.rs`, add to the struct (after `max_g_trace`):

```rust
    /// `max_t g` after each solve — non-increasing; read only by tests.
    #[allow(dead_code)]
    pub max_g_trace: Vec<f64>,
    /// Size of `T^est` solved at each iteration — non-monotone when a slack
    /// time is dropped then a violated maximum is re-added. Read only by tests.
    #[allow(dead_code)]
    pub active_set_trace: Vec<usize>,
```

- [ ] **Step 2: Populate it in the `refine` loop**

In `refine`, add `let mut active_set_trace = Vec::new();` next to `let mut max_g_trace = Vec::new();`. Inside the loop, after the empty-check and *before* assembling `rows`, record the size:

```rust
        active_set_trace.push(t_est.len());

        // Solve eq. 40 over the current candidate set T^est.
        let rows: Vec<ConicRows> = t_est
```

And add `active_set_trace,` to the `RefineOutcome { … }` returned on convergence (next to `max_g_trace,`).

- [ ] **Step 3: Verify the crate still builds and existing refine tests pass**

Run: `cargo test --lib algorithm::refine`
Expected: the existing 13 refine tests still PASS (the new field is additive).

- [ ] **Step 4: Write the real-`J2Roe` drop-then-readd test**

Add to the `#[cfg(test)] mod tests` in `src/algorithm/refine.rs`. The setup builds a **two-impulse reachable** target on the real `J2Roe` `Γ`, seeds `T^est` with a coarse set that includes a decoy slack time (to be dropped) but misses one optimal time (to be re-added), forcing ≥3 iterations with a non-monotone active set:

```rust
    #[test]
    fn refine_on_real_j2roe_runs_three_iters_with_drop_then_readd() {
        use crate::dynamics::{AbsoluteOrbit, J2Roe};

        // Real ill-conditioned worked-example dynamics (Γ δλ-row ~1e3, others ~1e-4).
        let chief = AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        let (t_i, t_f) = (0.0, 117_990.0);
        let dynamics = J2Roe::new(chief, t_i, t_f);
        // Coarse 300 s grid (394 points) — enough resolution, cheap to cache.
        let grid = TimeGrid::uniform(t_i, t_f, 300.0);
        let gammas: Vec<SMatrix<f64, N, M>> = grid.times().map(|t| dynamics.gamma(t)).collect();
        let cost = Piecewise::new(std::f64::consts::TAU / chief.mean_motion());

        // Two-impulse reachable target at two distinct grid times.
        let (ka, kb) = (120usize, 360usize); // t = 36000 s, 108000 s
        let ua = SVector::<f64, M>::new(0.6, -0.4, 0.5);
        let ub = SVector::<f64, M>::new(-0.3, 0.5, 0.4);
        let w = gammas[ka] * ua + gammas[kb] * ub;

        // Seed: one true time (ka), a decoy slack time (200), but NOT kb — so the
        // first solve drops the decoy and a later iteration re-adds a max near kb.
        let params = SolveParams::default();
        let out = refine(&cost, &grid, &gammas, &w, &params, vec![ka, 200], 50).unwrap();

        eprintln!(
            "[refine/j2roe] iters={} active_set_trace={:?} max_g_trace={:?}",
            out.iterations, out.active_set_trace, out.max_g_trace,
        );

        // Observably ≥3 iterations on the real ill-conditioned Γ.
        assert!(out.iterations >= 3, "iterations = {} (want ≥ 3)", out.iterations);
        // max_t g is non-increasing toward 1.
        for pair in out.max_g_trace.windows(2) {
            assert!(pair[1] <= pair[0] + 1e-6, "trace not non-increasing: {:?}", out.max_g_trace);
        }
        // Drop-then-readd: the active-set size dips below a predecessor and later grows again.
        let sizes = &out.active_set_trace;
        let dipped_then_grew = (1..sizes.len()).any(|i| {
            sizes[i] < sizes[i - 1] && sizes[i + 1..].iter().any(|&s| s > sizes[i])
        });
        assert!(dipped_then_grew, "expected a drop-then-readd in active set: {sizes:?}");
    }
```

- [ ] **Step 5: Run it and tune the seed via the printed traces if needed**

Run: `cargo test --lib refine_on_real_j2roe_runs_three_iters_with_drop_then_readd -- --nocapture`
Expected: PASS, with `iters ≥ 3` and an `active_set_trace` that dips then grows (e.g. `[2, 1, 2, …]` or `[2, 2, 1, 3, …]`).

**If `dipped_then_grew` is false** (the run converges without a visible re-add): inspect the printed `active_set_trace` and adjust the seed to force the pattern — options, in order: (i) add a second decoy slack time to the seed (e.g. `vec![ka, 150, 250]`) so a drop is guaranteed; (ii) move the seeded true time slightly off the optimum (e.g. `ka - 1`) so a neighboring max must be added; (iii) widen the grid resolution so the added local max lands on a fresh index. **Documented fallback** (only if a deterministic dip proves brittle across machines): drop the `dipped_then_grew` assertion and keep `iterations >= 3` + non-increasing `max_g_trace` on the real `J2Roe` `Γ` — this still satisfies the deferred item's core ("observably runs ≥3 iterations" on the real dynamics); leave a comment citing this plan step.

- [ ] **Step 6: Lint + commit**

Run: `cargo clippy --all-features --lib --tests -- -D warnings`
Expected: clean (the two `#[allow(dead_code)]` fields are read by this test, but the attribute is harmless and matches the `max_g_trace` precedent).

```bash
git add src/algorithm/refine.rs
git commit -m "test(phase5): refine on real J2Roe runs ≥3 iters with drop-then-readd (Phase-4 deferral)"
```

---

## Task 6: Deferred Phase-4 item (b) — `achieved > target` in the NotConverged test

Phase 4 deferred adding an `achieved > target` assertion to the iteration-cap test. Strengthen the existing `refine_reports_not_converged_at_iteration_cap`.

**Files:**
- Modify: `src/algorithm/refine.rs` (the `refine_reports_not_converged_at_iteration_cap` test)

- [ ] **Step 1: Bind and assert `achieved`/`target`**

In `src/algorithm/refine.rs`, change the `match err` arm of `refine_reports_not_converged_at_iteration_cap` from:

```rust
        match err {
            PlannerError::NotConverged { max_iters, .. } => assert_eq!(max_iters, 1),
            other => panic!("expected NotConverged, got {other:?}"),
        }
```

to:

```rust
        match err {
            PlannerError::NotConverged { max_iters, achieved, target } => {
                assert_eq!(max_iters, 1);
                // Not converged ⇒ the achieved max_t g must exceed the 1+eps_cost target.
                assert!(achieved > target, "achieved {achieved} should exceed target {target}");
            }
            other => panic!("expected NotConverged, got {other:?}"),
        }
```

- [ ] **Step 2: Run the test**

Run: `cargo test --lib refine_reports_not_converged_at_iteration_cap`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add src/algorithm/refine.rs
git commit -m "test(phase5): assert achieved > target in NotConverged cap test (Phase-4 deferral)"
```

---

## Task 7: Full gate, spec + memory update, PR

Run the complete CI gate locally, record the validated numbers in the spec, update memory, and open the PR.

**Files:**
- Modify: `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`
- Modify: `…/memory/spec-validation-status.md` and `…/memory/MEMORY.md`

- [ ] **Step 1: Run the full CI gate locally (must all pass)**

```bash
cargo fmt --all -- --check
cargo clippy --all-features --all-targets -- -D warnings
cargo build --all-features
cargo test --all-features
```
Expected: all green. Note the new test count (the prior total was 89 lib+integration tests after Phase 4; Phase 5 adds the worked-example file, the Hunter case, and the two refine tests). Capture the exact `test result: ok. N passed` lines.

- [ ] **Step 2: Mark Phase 5 done in the spec**

In `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`:
- Update the `Status:` line (§ top, line ~4) from "Phases 0–4 complete … Resume at Phase 5" to "Phases 0–5 complete … Resume at Phase 6 (Monte Carlo)."
- Replace the `### Phase 5 — Worked-example validation` body (line ~495) with a `✅ Done` entry recording: the files added (`tests/worked_example.rs`, `examples/mdot.rs` filled in); the **confirmed `a_c` scaling** (`w_nd = w_metres / a_c`, Δv in m/s, λ scaled by `a_c` ⇒ direction-only assertion); the **observed** maneuvers/total_dv/iterations/residual vs §7 (paste Task 1's numbers and any band that needed a documented tolerance widening); the Fig. 7 `max_t g` value; the Hunter cross-check result and which `u₀` reading reproduced it; and the two Phase-4 deferrals closed (the real-`J2Roe` drop-then-readd test, the `achieved > target` assertion). Keep the prose style of the Phase 1–4 entries.

- [ ] **Step 3: Update memory**

In `…/memory/spec-validation-status.md`, append a Phase-5 clause to the status (Phase 5 ✅: worked example reproduces Table IV within §7 bands; `examples/mdot.rs` + `tests/worked_example.rs` landed; next = Phase 6 Monte Carlo). Update the matching one-line hook in `…/memory/MEMORY.md`. Convert any relative dates to absolute (today = 2026-06-18).

- [ ] **Step 4: Commit the docs**

```bash
git add docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md
git commit -m "docs(spec): mark Phase 5 (worked-example validation) complete; resume at Phase 6"
```
(Memory files live outside the repo working tree — they are saved by the Write tool, not committed.)

- [ ] **Step 5: Push and open the PR**

```bash
git push -u origin phase5-worked-example
gh pr create --base main --title "Phase 5 — worked-example validation (Table III -> Table IV, Fig. 7)" \
  --body "Reproduces the §7 worked example end-to-end (3 maneuvers, ≈82.4 mm/s, ~3 iters, residual < 0.01%) as a CI-enforced integration test; adds examples/mdot.rs (Table IV report + Fig. 7 contact curve); closes the two Phase-4-deferred refine tests. Bands per spec Risk 3/8."
```
Then use the `superpowers:finishing-a-development-branch` skill to confirm CI is green and merge.

---

## Self-Review

**1. Spec coverage (§ Phase 5 + §7):**
- "Encode Table III + params" → Task 1 (`worked_example_inputs`, `SolveParams::default`). ✓
- "reproduce the published result" / "all §7 worked-example assertions pass within stated bands" → Task 2 (3 maneuvers, times, per-component Δv, total_dv, iterations, residual, λ direction). ✓
- Fig. 7 contact-function curve → Task 3 (computed over the grid, CSV behind `validation`, `max_t g ≤ 1.01` self-check). ✓
- Second worked example (Hunter, §7) → Task 4. ✓ (marked secondary, as the spec directs).
- Phase-4 deferrals (real-`J2Roe` ≥3-iter drop-then-readd; `achieved > target`) → Tasks 5–6. ✓
- `a_c` scaling (§5.5) → Global Constraints + Task 1 (`/ A_C`, mm/s reporting, λ direction). ✓
- CI gate green (exit criterion) → Task 7 Step 1. ✓

**2. Placeholder scan:** No "TBD/handle edge cases/similar to Task N". Every code step shows complete code; every run step gives the exact command + expected result. The two characterization decision points (Task 1 Step 4, Task 4 Step 2, Task 5 Step 5) are explicit triage procedures, not placeholders. ✓

**3. Type consistency:** `worked_example_inputs()` returns `(J2Roe, Piecewise, Pseudostate, TimeGrid, SolveParams)` and is consumed verbatim in Tasks 1–2. `A_C`/`W_METRES` are reused consistently (and deliberately duplicated into `examples/mdot.rs` per the documented Rust constraint). `Solution` field names (`maneuvers`, `total_dv`, `iterations`, `residual`, `lambda`) match `src/types.rs`. `RefineOutcome` gains `active_set_trace` consistently across the struct, the loop, and the return (Task 5). `solve`/`AbsoluteOrbit::new`/`J2Roe::new`/`Piecewise::new`/`TimeGrid::uniform` signatures match the on-disk verbatim. `dv` reported as `* 1e3` for mm/s throughout. ✓

**Risk note carried forward:** the dominant uncertainty is whether clarabel reproduces §7 to within the chosen bands on the ill-conditioned `Γ`. The characterization-first structure surfaces any miss as data before any band is committed; structural misses route to systematic-debugging, tolerance misses to a documented band widening — never the reverse.
