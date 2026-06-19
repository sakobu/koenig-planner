# Koenig Planner — Phase 6: Monte Carlo Harness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Monte Carlo validation harness (`src/bin/monte_carlo.rs`) reproducing Fig. 8 (Algorithm-2 iteration distributions for 2/6/10-time inits) and Fig. 9 (compute time vs discretization `|T|`, 10→10⁶) for the *proposed* algorithm, plus a seeded CI test asserting paper-independent invariants.

**Architecture:** One feature-gated binary owns sampling, both sweeps, CSV/PNG output, and bin-internal unit tests for its pure helpers. A separate seeded integration test (`tests/monte_carlo.rs`) asserts robust invariants via the public `solve` API. No harness code is added to the library crate; the sampling *convention* is documented so the bin and the test agree, but they sample independently (the asserted invariants hold for any seed).

**Tech Stack:** Rust, `nalgebra` (`SVector<f64,6>`), `clarabel` (already wired through `solve`), `rand 0.10` + `rand_distr 0.6` (seeded portable `StdRng` + Gaussian, NEW, optional), `csv 1.4` + `plotters 0.3` (already optional). All new deps live behind the existing `validation` feature.

## Global Constraints

- **CI gate (every task ends green):** `cargo fmt --check`, `cargo clippy --all-features -- -D warnings`, `cargo build --all-features`, `cargo test --all-features` — all `--all-features`. The Linux runner installs `libfontconfig1-dev` for `plotters`.
- **Formatting:** the code snippets below are written for readability, **not** pre-wrapped to rustfmt's 100-col output (several struct literals and `assert!` macros exceed one line once pasted). Always run `cargo fmt --all` after pasting a snippet — the gate commands below begin with it (`cargo fmt --all && cargo fmt --check && …`), which formats and then verifies CI-cleanliness. CI itself runs only `cargo fmt --check`, so committed code must already be formatted.
- **Edition 2021, `rust-version = "1.92"`.** Existing deps: `nalgebra 0.35`, `clarabel 0.11`, `thiserror 2.0`, `approx 0.5` (dev), `csv 1.4` + `plotters 0.3` (optional). New: `rand 0.10`, `rand_distr 0.6` (optional).
- **Determinism:** use `StdRng` (ChaCha-based, portable macOS ↔ Linux CI) seeded from a documented constant — **never** `SmallRng` (non-portable) or `rand::rng()`/thread RNG.
- **Scaling convention (spec §5.5):** `a_c = 25_000e3` m; feed `solve` a dimensionless `w_nd = w_metres / a_c`; Δv comes out in m/s.
- **Sampling convention (spec §6 Phase 6):** each of the 6 ROE components `~ Normal(0, σ = 1000 m)` (metre-scaled `w`), then `÷ a_c`.
- **Fixed MC problem = the worked example (Table III):** chief `AbsoluteOrbit::new(a_c, 0.7, 40°, 358°, 0°, 180°)` (angles in radians); `t_i = 0`, `t_f = 117_990` s; cost `Piecewise::new(TAU / chief.mean_motion())` (eq. 49); `SolveParams::default()` except `n_init`.
- **Public API (do not change signatures):** `solve(&dyn, &cost, w: Pseudostate, grid: TimeGrid, &params) -> Result<Solution, PlannerError>`; `Solution { maneuvers, total_dv, iterations, residual, lambda }`; `TimeGrid::uniform(t_i, t_f, dt)` + `.len()`; `SolveParams { n_coarse, n_init, eps_cost, eps_remove, q }` + `::default()`; `J2Roe::new(chief, t_i, t_f)`; `AbsoluteOrbit::new(a, e, i, Ω, ω, M)` + `.mean_motion()`; `Piecewise::new(period_s)`; `Pseudostate = SVector<f64,6>` with `::from_row_slice(&[f64;6])`.
- **Validation stance:** report our distributions and compare to the paper's `4.90 / 3.99 / 3.31` means as *reference*, NOT pass/fail. The CI test asserts only paper-independent invariants. Per the Phase-4/5 band methodology: implement, run to observe, then lock bands with margin (bands, not bit-equality — risk #3).

## File Structure

| File | Responsibility |
|---|---|
| `Cargo.toml` (modify) | Add `rand`/`rand_distr` optional deps; extend `validation = [...]`. |
| `src/bin/monte_carlo.rs` (rewrite the stub) | The whole harness: constants, fixed-problem builders, seeded sampler, Fig. 8 sweep + summary + CSV, Fig. 9 sweep + CSV, optional PNGs, gated `main`, and bin-internal `#[cfg(test)]` unit tests for the pure helpers. Everything real is under `#[cfg(feature = "validation")]`. |
| `tests/monte_carlo.rs` (create) | Seeded integration test asserting MC invariants via the public `solve` API. Gated `#![cfg(feature = "validation")]`. |

**Gating shape (why):** a `src/bin/*.rs` target always compiles under a plain `cargo build`, so it cannot reference `rand`/`csv`/`plotters` (only present under `validation`) at file scope. The entire real harness therefore lives inside `#[cfg(feature = "validation")] mod harness { ... }`; a `#[cfg(not(feature = "validation"))] fn main()` prints a one-line hint. This mirrors how `examples/mdot.rs` gates its CSV block. Bin-internal `#[cfg(test)]` unit tests run under `cargo test --features validation` (CI uses `--all-features`).

---

### Task 1: Dependencies, feature wiring, and the gated bin skeleton

**Files:**
- Modify: `Cargo.toml`
- Rewrite: `src/bin/monte_carlo.rs` (currently a 5-line stub)

**Interfaces:**
- Produces (for later tasks, all inside `mod harness`): `const A_C: f64`, `const SIGMA_M: f64`, `const SEED: u64`, `const T_I: f64`, `const T_F: f64`, `const GRID_DT: f64`; `fn worked_example_chief() -> AbsoluteOrbit`, `fn worked_example_dynamics() -> J2Roe`, `fn worked_example_cost() -> Piecewise`; `fn main()` (delegates to the gated harness).

- [ ] **Step 1: Add the optional deps and extend the feature in `Cargo.toml`**

Add two dependency tables next to the existing `[dependencies.csv]` / `[dependencies.plotters]` tables:

```toml
[dependencies.rand]
version = "0.10"
optional = true

[dependencies.rand_distr]
version = "0.6"
optional = true
```

And replace the `[features]` section:

```toml
[features]
validation = ["dep:csv", "dep:plotters", "dep:rand", "dep:rand_distr"]
```

- [ ] **Step 2: Verify the deps resolve to the pinned versions**

Run: `cargo update -p rand --precise 0.10.1 && cargo update -p rand_distr --precise 0.6.0 && cargo tree -e features -i rand 2>/dev/null | head -5`
Expected: resolves `rand v0.10.1` and `rand_distr v0.6.0` (both in the local cargo cache, so this works offline). If `cargo update --precise` complains the package isn't a dependency yet, run a no-op `cargo build --features validation` first to populate `Cargo.lock`, then re-run.

- [ ] **Step 3: Rewrite `src/bin/monte_carlo.rs` to the gated skeleton**

Replace the entire file with:

```rust
//! Monte Carlo harness — Fig. 8 (Algorithm-2 iteration distributions) and Fig. 9
//! (compute time vs discretization |T|) for the *proposed* algorithm on the
//! worked-example problem (Table III chief, eq. 49 cost).
//!
//! Build & run with the `validation` feature (needs rand/csv/plotters):
//!
//!   cargo run --features validation --bin monte_carlo            # both sweeps
//!   cargo run --features validation --bin monte_carlo -- fig8    # Fig. 8 only
//!   cargo run --features validation --bin monte_carlo -- fig9    # Fig. 9 only
//!
//! Validation stance (spec §6 Phase 6): we REPORT our iteration distributions and
//! compare to the paper's 4.90/3.99/3.31 means as *reference*, not bit-reproduce the
//! paper's (internally inconsistent) figures. See tests/monte_carlo.rs for the
//! asserted, paper-independent invariants.

#[cfg(not(feature = "validation"))]
fn main() {
    eprintln!("monte_carlo: rebuild with `--features validation` (needs rand, csv, plotters).");
}

#[cfg(feature = "validation")]
fn main() {
    harness::main();
}

#[cfg(feature = "validation")]
mod harness {
    use koenig_planner::cost::Piecewise;
    use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
    use std::f64::consts::TAU;

    /// Chief semimajor axis a_c [m] — the I/O scaling factor (spec §5.5).
    pub const A_C: f64 = 25_000e3;
    /// Per-ROE Gaussian std, metre-scaled (σ = 1 km; spec §6 Phase 6).
    pub const SIGMA_M: f64 = 1000.0;
    /// Documented constant seed (portable StdRng) — "koenig" in hex-ish.
    pub const SEED: u64 = 0x6F_656E_6967;
    /// Worked-example window [s].
    pub const T_I: f64 = 0.0;
    pub const T_F: f64 = 117_990.0;
    /// Fig. 8 grid step [s] (Table III 30 s grid → 3934 candidate times).
    pub const GRID_DT: f64 = 30.0;

    /// Table III chief mean absolute orbit (angles in radians).
    pub fn worked_example_chief() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            A_C,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        )
    }

    /// J2 mean-ROE dynamics for the worked-example window.
    pub fn worked_example_dynamics() -> J2Roe {
        J2Roe::new(worked_example_chief(), T_I, T_F)
    }

    /// eq. 49 piecewise cost (FaceMax in 2-hr perigee windows, Norm2 elsewhere).
    pub fn worked_example_cost() -> Piecewise {
        Piecewise::new(TAU / worked_example_chief().mean_motion())
    }

    pub fn main() {
        let dynamics = worked_example_dynamics();
        let cost = worked_example_cost();
        let _ = (&dynamics, &cost); // wired into the sweeps in later tasks
        println!("koenig-planner Monte Carlo harness (Phase 6)");
        println!("  problem            : worked example (Table III chief, eq. 49 cost)");
        println!("  window [s]         : [{T_I}, {T_F}]");
        println!("  mean motion [rad/s]: {:.6e}", worked_example_chief().mean_motion());
        println!("  Fig. 8 grid        : dt = {GRID_DT} s");
        println!("  seed               : {SEED:#x}");
    }
}
```

- [ ] **Step 4: Verify the default build still works (stub main)**

Run: `cargo run --bin monte_carlo 2>&1`
Expected: prints `monte_carlo: rebuild with --features validation ...` (the `#[cfg(not(...))]` path; goes to stderr). `cargo build` (no features) compiles clean.

- [ ] **Step 5: Verify the validation build runs the real main**

Run: `cargo run --features validation --bin monte_carlo`
Expected: prints the problem header with a nonzero mean motion (≈ `1.597e-4` rad/s) and the seed.

- [ ] **Step 6: Run the full gate**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green (existing 105 tests still pass; no new tests yet).

- [ ] **Step 7: Commit**

```bash
git add Cargo.toml Cargo.lock src/bin/monte_carlo.rs
git commit -m "feat(mc): add rand/rand_distr deps + gated monte_carlo skeleton"
```

---

### Task 2: Seeded Gaussian pseudostate sampler (TDD)

**Files:**
- Modify: `src/bin/monte_carlo.rs` (inside `mod harness`)

**Interfaces:**
- Consumes: `A_C`, `SIGMA_M` (Task 1).
- Produces: `fn sample_pseudostates(n: usize, seed: u64) -> Vec<Pseudostate>` — `n` dimensionless `w_nd` targets.

- [ ] **Step 1: Add the imports and a failing test**

At the top of `mod harness`, extend the `use` block:

```rust
    use koenig_planner::Pseudostate;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal};
```

Add a bin-internal test module at the bottom of `mod harness` (this only compiles under `--features validation`, so no extra gate is needed beyond the module's own `#[cfg(test)]`):

```rust
    #[cfg(test)]
    mod tests {
        use super::*;
        use approx::assert_abs_diff_eq;

        #[test]
        fn sampler_is_deterministic_and_well_scaled() {
            let a = sample_pseudostates(200, SEED);
            let b = sample_pseudostates(200, SEED);
            assert_eq!(a.len(), 200);
            // Determinism: same seed -> identical samples.
            for (x, y) in a.iter().zip(&b) {
                assert_eq!(x, y);
            }
            // Convention: components ~ Normal(0, SIGMA_M / A_C); never near-zero norm.
            let expected_sd = SIGMA_M / A_C;
            let flat: Vec<f64> = a.iter().flat_map(|w| w.iter().copied()).collect();
            let mean = flat.iter().sum::<f64>() / flat.len() as f64;
            let var = flat.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / flat.len() as f64;
            assert_abs_diff_eq!(mean, 0.0, epsilon = expected_sd * 0.15);
            assert_abs_diff_eq!(var.sqrt(), expected_sd, epsilon = expected_sd * 0.15);
            assert!(a.iter().all(|w| w.norm() > 0.0));
        }

        #[test]
        fn different_seeds_differ() {
            let a = sample_pseudostates(8, SEED);
            let b = sample_pseudostates(8, SEED + 1);
            assert!(a.iter().zip(&b).any(|(x, y)| x != y));
        }
    }
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --features validation --bin monte_carlo sampler -- --nocapture`
Expected: FAIL — `cannot find function sample_pseudostates`.

- [ ] **Step 3: Implement the sampler**

Add to `mod harness` (after the builders):

```rust
    /// `n` random target pseudostates as dimensionless `w_nd`: each of the 6 ROE
    /// components `~ Normal(0, σ = SIGMA_M metres)`, then divided by `a_c`
    /// (spec §6 Phase 6 sampling convention). `StdRng` is portable, so a fixed
    /// `seed` yields identical samples on every platform.
    pub fn sample_pseudostates(n: usize, seed: u64) -> Vec<Pseudostate> {
        let mut rng = StdRng::seed_from_u64(seed);
        let normal = Normal::new(0.0_f64, SIGMA_M).expect("σ > 0 is a valid normal");
        (0..n)
            .map(|_| {
                let mut comp = [0.0_f64; 6];
                for c in comp.iter_mut() {
                    *c = normal.sample(&mut rng);
                }
                Pseudostate::from_row_slice(&comp) / A_C
            })
            .collect()
    }
```

`sample_pseudostates` is not yet called by `main`, so add a transient `#[allow(dead_code)]` directly above `pub fn sample_pseudostates` (removed in Task 3 when the Fig. 8 sweep calls it — the established Phase-3/4 transient-allow pattern; `#[expect]` mis-fires under `cfg(test)`).

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test --features validation --bin monte_carlo sampler different_seeds -- --nocapture`
Expected: PASS (both `sampler_is_deterministic_and_well_scaled` and `different_seeds_differ`).

- [ ] **Step 5: Run the full gate**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green. (The `#[allow(dead_code)]` keeps clippy quiet about the not-yet-wired sampler.)

- [ ] **Step 6: Commit**

```bash
git add src/bin/monte_carlo.rs
git commit -m "feat(mc): seeded portable Gaussian pseudostate sampler (TDD)"
```

---

### Task 3: Fig. 8 sweep — iteration distributions, summary, CSV

**Files:**
- Modify: `src/bin/monte_carlo.rs` (inside `mod harness`)

**Interfaces:**
- Consumes: `sample_pseudostates` (Task 2), `worked_example_dynamics`/`worked_example_cost` (Task 1), `T_I`/`T_F`/`GRID_DT`/`SEED`.
- Produces: `struct Fig8Row { n_init, sample, iterations, residual, total_dv }`; `struct Fig8Stat { n_init, n, mean_iters, frac_le8, max_iters, max_residual }`; `fn run_fig8(&D, &C, &[Pseudostate], &[usize]) -> (Vec<Fig8Row>, usize)`; `fn summarize_fig8(&[Fig8Row], &[usize]) -> Vec<Fig8Stat>`; `fn write_fig8_csv(&str, &[Fig8Row]) -> csv::Result<()>`; `fn fig8(&D, &C)`; `const N_MC: usize`, `const N_INITS: [usize; 3]`, `const PAPER_MEANS: [f64; 3]`.

- [ ] **Step 1: Extend the imports and add a failing smoke test**

Extend the `use` block in `mod harness`:

```rust
    use koenig_planner::{solve, CostModel, Dynamics, SolveParams, TimeGrid};
```

Add to the `#[cfg(test)] mod tests` block:

```rust
        #[test]
        fn fig8_sweep_produces_paired_rows_and_low_residual() {
            let dynamics = worked_example_dynamics();
            let cost = worked_example_cost();
            let ws = sample_pseudostates(3, SEED);
            let n_inits = [2usize, 6];
            let (rows, failures) = run_fig8(&dynamics, &cost, &ws, &n_inits);
            assert_eq!(failures, 0, "no solve should fail on the worked-example problem");
            assert_eq!(rows.len(), 3 * 2, "one row per (n_init, sample)");
            for r in &rows {
                assert!(r.residual < 1e-3, "row residual {:.3e} too high", r.residual);
                assert!((1..=50).contains(&r.iterations));
            }
            let stats = summarize_fig8(&rows, &n_inits);
            assert_eq!(stats.len(), 2);
            assert!(stats.iter().all(|s| s.n == 3));
        }
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --features validation --bin monte_carlo fig8_sweep -- --nocapture`
Expected: FAIL — `cannot find function run_fig8` / `summarize_fig8`.

- [ ] **Step 3: Implement the sweep, summary, and CSV writer**

Add to `mod harness`. First the constants near the top:

```rust
    /// Fig. 8 sample count (paper: 200).
    pub const N_MC: usize = 200;
    /// Initialization candidate counts swept by Fig. 8.
    pub const N_INITS: [usize; 3] = [2, 6, 10];
    /// Paper's reported mean iteration counts for N_INITS (reference, not a target).
    pub const PAPER_MEANS: [f64; 3] = [4.90, 3.99, 3.31];
```

Then the types and functions:

```rust
    /// One Fig. 8 sample outcome.
    #[derive(Clone, Copy)]
    pub struct Fig8Row {
        pub n_init: usize,
        pub sample: usize,
        pub iterations: usize,
        pub residual: f64,
        pub total_dv: f64,
    }

    /// Per-`n_init` summary statistics.
    pub struct Fig8Stat {
        pub n_init: usize,
        pub n: usize,
        pub mean_iters: f64,
        pub frac_le8: f64,
        pub max_iters: usize,
        pub max_residual: f64,
    }

    /// Run `solve` for every `(n_init, w)` pair on the fixed 30 s grid; collect
    /// per-sample outcomes and a count of solver failures (Phase 5b ⇒ expect 0).
    pub fn run_fig8<D: Dynamics, C: CostModel>(
        dynamics: &D,
        cost: &C,
        ws: &[Pseudostate],
        n_inits: &[usize],
    ) -> (Vec<Fig8Row>, usize) {
        let grid = TimeGrid::uniform(T_I, T_F, GRID_DT);
        let mut rows = Vec::with_capacity(ws.len() * n_inits.len());
        let mut failures = 0usize;
        for &n_init in n_inits {
            let params = SolveParams { n_init, ..SolveParams::default() };
            for (sample, &w) in ws.iter().enumerate() {
                match solve(dynamics, cost, w, grid, &params) {
                    Ok(sol) => rows.push(Fig8Row {
                        n_init,
                        sample,
                        iterations: sol.iterations,
                        residual: sol.residual,
                        total_dv: sol.total_dv,
                    }),
                    Err(_) => failures += 1,
                }
            }
        }
        (rows, failures)
    }

    /// Group rows by `n_init` and compute mean iterations, fraction ≤ 8 iterations,
    /// max iterations, and max residual.
    pub fn summarize_fig8(rows: &[Fig8Row], n_inits: &[usize]) -> Vec<Fig8Stat> {
        n_inits
            .iter()
            .map(|&n_init| {
                let group: Vec<&Fig8Row> = rows.iter().filter(|r| r.n_init == n_init).collect();
                let n = group.len();
                let denom = n.max(1) as f64;
                let mean_iters = group.iter().map(|r| r.iterations as f64).sum::<f64>() / denom;
                let frac_le8 = group.iter().filter(|r| r.iterations <= 8).count() as f64 / denom;
                let max_iters = group.iter().map(|r| r.iterations).max().unwrap_or(0);
                let max_residual = group.iter().map(|r| r.residual).fold(0.0, f64::max);
                Fig8Stat { n_init, n, mean_iters, frac_le8, max_iters, max_residual }
            })
            .collect()
    }

    /// Write the per-sample Fig. 8 rows to `path` as CSV.
    pub fn write_fig8_csv(path: &str, rows: &[Fig8Row]) -> csv::Result<()> {
        let mut w = csv::Writer::from_path(path)?;
        w.write_record(["n_init", "sample", "iterations", "residual", "total_dv"])?;
        for r in rows {
            w.write_record(&[
                r.n_init.to_string(),
                r.sample.to_string(),
                r.iterations.to_string(),
                format!("{:.6e}", r.residual),
                format!("{:.9e}", r.total_dv),
            ])?;
        }
        w.flush()?;
        Ok(())
    }

    /// Fig. 8 driver: sample, sweep, summarize, print, and write the CSV.
    pub fn fig8<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C) {
        let ws = sample_pseudostates(N_MC, SEED);
        let (rows, failures) = run_fig8(dynamics, cost, &ws, &N_INITS);
        let stats = summarize_fig8(&rows, &N_INITS);

        println!("\nFig. 8 — Algorithm-2 iteration distribution ({N_MC} samples/init)");
        println!(
            "  {:>6}  {:>5}  {:>10}  {:>8}  {:>11}  {:>12}",
            "n_init", "n", "mean_iters", "frac<=8", "max_iters", "max_residual"
        );
        for (s, paper) in stats.iter().zip(PAPER_MEANS.iter()) {
            println!(
                "  {:>6}  {:>5}  {:>10.3}  {:>8.3}  {:>11}  {:>12.2e}   (paper {:.2})",
                s.n_init, s.n, s.mean_iters, s.frac_le8, s.max_iters, s.max_residual, paper
            );
        }
        if failures > 0 {
            println!("  WARNING: {failures} solve(s) failed (expected 0).");
        }

        let path = "target/fig8_iterations.csv";
        match write_fig8_csv(path, &rows) {
            Ok(()) => println!("  rows written         : {path} ({} rows)", rows.len()),
            Err(e) => eprintln!("  CSV write failed     : {e}"),
        }
    }
```

Remove the transient `#[allow(dead_code)]` from `sample_pseudostates` (now called by `fig8`).

- [ ] **Step 4: Wire `fig8` into `main` and remove the placeholder line**

In `mod harness::main`, replace the `let _ = (&dynamics, &cost);` placeholder line and the four trailing header `println!`s with arg handling that calls `fig8` (Fig. 9 is added in Task 4):

```rust
    pub fn main() {
        let which = std::env::args().nth(1);
        if let Some(arg) = which.as_deref() {
            if arg != "fig8" && arg != "fig9" {
                eprintln!("usage: monte_carlo [fig8|fig9]   (default: both)");
                std::process::exit(2);
            }
        }
        std::fs::create_dir_all("target").ok();
        let dynamics = worked_example_dynamics();
        let cost = worked_example_cost();

        println!("koenig-planner Monte Carlo harness (Phase 6)  seed={SEED:#x}");
        let run_8 = matches!(which.as_deref(), None | Some("fig8"));
        if run_8 {
            fig8(&dynamics, &cost);
        }
    }
```

- [ ] **Step 5: Run the smoke test and the Fig. 8 run**

Run: `cargo test --features validation --bin monte_carlo fig8_sweep -- --nocapture`
Expected: PASS.

Run: `cargo run --features validation --bin monte_carlo -- fig8`
Expected: prints the 3-row summary table (one row per `n_init` with mean iterations, `frac<=8 ≈ 1.000`, tiny `max_residual`), writes `target/fig8_iterations.csv` with `200*3 = 600` data rows. Observe and note the three mean iteration values (compared to 4.90/3.99/3.31).

- [ ] **Step 6: Run the full gate**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green.

- [ ] **Step 7: Commit**

```bash
git add src/bin/monte_carlo.rs
git commit -m "feat(mc): Fig. 8 iteration-distribution sweep + summary + CSV"
```

---

### Task 4: Fig. 9 sweep — compute time vs |T|, CSV

**Files:**
- Modify: `src/bin/monte_carlo.rs` (inside `mod harness`)

**Interfaces:**
- Consumes: `solve`, `TimeGrid`, `SolveParams`, `Dynamics`, `CostModel`, `Pseudostate`, `T_I`/`T_F`, `sample_pseudostates`.
- Produces: `struct Fig9Row { grid_len, dt, seconds, iterations, residual }`; `fn run_fig9(&D, &C, Pseudostate, &[usize]) -> Vec<Fig9Row>`; `fn write_fig9_csv(&str, &[Fig9Row]) -> csv::Result<()>`; `fn fig9(&D, &C)`; `const FIG9_SIZES: [usize; 6]`.

- [ ] **Step 1: Add `std::time::Instant` and a failing smoke test**

Extend the `use` block in `mod harness`:

```rust
    use std::time::Instant;
```

Add to `#[cfg(test)] mod tests`:

```rust
        #[test]
        fn fig9_sweep_times_each_size() {
            let dynamics = worked_example_dynamics();
            let cost = worked_example_cost();
            let w = sample_pseudostates(1, SEED)[0];
            let sizes = [10usize, 100];
            let rows = run_fig9(&dynamics, &cost, w, &sizes);
            assert_eq!(rows.len(), 2);
            assert!(rows.iter().all(|r| r.seconds >= 0.0 && r.grid_len >= 2));
            // Finer grid is at least as large in point count.
            assert!(rows[1].grid_len >= rows[0].grid_len);
        }
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --features validation --bin monte_carlo fig9_sweep -- --nocapture`
Expected: FAIL — `cannot find function run_fig9`.

- [ ] **Step 3: Implement the Fig. 9 sweep and CSV writer**

Add to `mod harness`. The size schedule near the constants:

```rust
    /// Fig. 9 grid sizes (10 → 10⁶). 10⁶ is ~150 MB Γ cache / multi-second; documented.
    pub const FIG9_SIZES: [usize; 6] = [10, 100, 1_000, 10_000, 100_000, 1_000_000];
```

Then:

```rust
    /// One Fig. 9 timing outcome.
    #[derive(Clone, Copy)]
    pub struct Fig9Row {
        pub grid_len: usize,
        pub dt: f64,
        pub seconds: f64,
        pub iterations: usize,
        pub residual: f64,
    }

    /// Time `solve` once per requested grid size on the fixed window. For each size:
    /// `dt = (t_f - t_i)/(n-1)`, one warmup solve (discarded), one timed solve. The
    /// actual `grid.len()` is recorded (it may differ from `n` by ±1 due to rounding
    /// in `TimeGrid::len`). Timing shape is `w`-independent; use a single fixed `w`.
    pub fn run_fig9<D: Dynamics, C: CostModel>(
        dynamics: &D,
        cost: &C,
        w: Pseudostate,
        sizes: &[usize],
    ) -> Vec<Fig9Row> {
        let params = SolveParams::default();
        sizes
            .iter()
            .map(|&n| {
                let dt = (T_F - T_I) / (n.max(2) - 1) as f64;
                let grid = TimeGrid::uniform(T_I, T_F, dt);
                let _ = solve(dynamics, cost, w, grid, &params); // warmup
                let start = Instant::now();
                let result = solve(dynamics, cost, w, grid, &params);
                let seconds = start.elapsed().as_secs_f64();
                match result {
                    Ok(s) => Fig9Row {
                        grid_len: grid.len(),
                        dt,
                        seconds,
                        iterations: s.iterations,
                        residual: s.residual,
                    },
                    Err(_) => Fig9Row {
                        grid_len: grid.len(),
                        dt,
                        seconds,
                        iterations: 0,
                        residual: f64::NAN,
                    },
                }
            })
            .collect()
    }

    /// Write the Fig. 9 timing rows to `path` as CSV.
    pub fn write_fig9_csv(path: &str, rows: &[Fig9Row]) -> csv::Result<()> {
        let mut w = csv::Writer::from_path(path)?;
        w.write_record(["grid_len", "dt_s", "seconds", "iterations", "residual"])?;
        for r in rows {
            w.write_record(&[
                r.grid_len.to_string(),
                format!("{:.6e}", r.dt),
                format!("{:.6e}", r.seconds),
                r.iterations.to_string(),
                format!("{:.3e}", r.residual),
            ])?;
        }
        w.flush()?;
        Ok(())
    }

    /// Fig. 9 driver: time the Table III `w` across `FIG9_SIZES`, print, write CSV.
    pub fn fig9<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C) {
        let w = sample_pseudostates(1, SEED)[0];
        println!("\nFig. 9 — solve time vs |T| (10⁶ is multi-second / ~150 MB)");
        let rows = run_fig9(dynamics, cost, w, &FIG9_SIZES);
        println!(
            "  {:>10}  {:>12}  {:>10}  {:>6}  {:>10}",
            "grid_len", "dt_s", "seconds", "iters", "residual"
        );
        for r in &rows {
            println!(
                "  {:>10}  {:>12.4e}  {:>10.4}  {:>6}  {:>10.2e}",
                r.grid_len, r.dt, r.seconds, r.iterations, r.residual
            );
        }
        let path = "target/fig9_timing.csv";
        match write_fig9_csv(path, &rows) {
            Ok(()) => println!("  rows written         : {path} ({} rows)", rows.len()),
            Err(e) => eprintln!("  CSV write failed     : {e}"),
        }
    }
```

- [ ] **Step 4: Wire `fig9` into `main`**

In `mod harness::main`, after the `if run_8 { fig8(...) }` block add:

```rust
        let run_9 = matches!(which.as_deref(), None | Some("fig9"));
        if run_9 {
            fig9(&dynamics, &cost);
        }
```

- [ ] **Step 5: Run the smoke test and a quick Fig. 9 run**

Run: `cargo test --features validation --bin monte_carlo fig9_sweep -- --nocapture`
Expected: PASS.

Run: `cargo run --release --features validation --bin monte_carlo -- fig9`
Expected: a 6-row table; `seconds` ≈ flat for `grid_len ≤ 10⁴` then rising, with the `10⁶` row taking the longest (seconds) and the largest memory. All `residual` tiny except possibly the smallest grids (few candidate times can't reconstruct `w` — that is expected for `|T|=10` and does not affect the timing measurement). Note the constant-then-linear shape. (Use `--release`; the 10⁶ point is slow in debug.)

- [ ] **Step 6: Run the full gate**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green. (The Fig. 9 *test* uses only sizes 10/100, so `cargo test` stays fast — the 10⁶ point runs only via `cargo run`.)

- [ ] **Step 7: Commit**

```bash
git add src/bin/monte_carlo.rs
git commit -m "feat(mc): Fig. 9 compute-time-vs-|T| sweep + CSV"
```

---

### Task 5: Optional PNGs — iteration CDF + log-log timing (plotters)

**Files:**
- Modify: `src/bin/monte_carlo.rs` (inside `mod harness`)

**Interfaces:**
- Consumes: `Fig8Row`/`Fig8Stat` not needed; raw `iterations` per `n_init`; `Fig9Row`.
- Produces: `fn empirical_cdf(&[usize]) -> Vec<(f64, f64)>`; `fn plot_fig8_cdf(&str, &[(usize, Vec<(f64, f64)>)]) -> Result<(), Box<dyn std::error::Error>>`; `fn plot_fig9_timing(&str, &[Fig9Row]) -> Result<(), Box<dyn std::error::Error>>`. Wires PNG writes into `fig8`/`fig9`.

> **plotters is the one external-API risk in this plan.** The code below is written for the resolved `plotters 0.3.7` (default features: `bitmap_backend` + `bitmap_encoder` give PNG output; CI installs `libfontconfig1-dev` for text). If a signature drifts, adjust against the cached `plotters-0.3.7` source — the CSVs from Tasks 3–4 are the actual deliverable; the PNGs are best-effort per spec ("optional `plotters` PNGs").

- [ ] **Step 1: Add a failing test for the CDF helper**

Add to `#[cfg(test)] mod tests`:

```rust
        #[test]
        fn empirical_cdf_is_monotone_and_ends_at_one() {
            let pts = empirical_cdf(&[3, 3, 4, 5]);
            // Left anchor at 0, then steps at the distinct values, ending at 1.0.
            assert_eq!(pts.first().unwrap().1, 0.0);
            assert!((pts.last().unwrap().1 - 1.0).abs() < 1e-12);
            // (value, fraction <= value): (2,0), (3,0.5), (4,0.75), (5,1.0).
            assert_eq!(pts, vec![(2.0, 0.0), (3.0, 0.5), (4.0, 0.75), (5.0, 1.0)]);
            // Monotone non-decreasing in the fraction.
            assert!(pts.windows(2).all(|w| w[1].1 >= w[0].1));
        }
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --features validation --bin monte_carlo empirical_cdf -- --nocapture`
Expected: FAIL — `cannot find function empirical_cdf`.

- [ ] **Step 3: Implement the CDF helper and the two plot functions**

Add `use plotters::prelude::*;` to the `mod harness` `use` block. Then add:

```rust
    /// Empirical CDF of iteration counts as `(value, fraction ≤ value)` over the
    /// distinct sorted values, anchored at `(min-1, 0)` so the curve starts at 0.
    pub fn empirical_cdf(counts: &[usize]) -> Vec<(f64, f64)> {
        assert!(!counts.is_empty(), "empirical_cdf needs at least one sample");
        let n = counts.len() as f64;
        let mut sorted = counts.to_vec();
        sorted.sort_unstable();
        let mut pts = vec![((*sorted.first().unwrap() as f64) - 1.0, 0.0)];
        let mut i = 0;
        while i < sorted.len() {
            let v = sorted[i];
            let mut j = i;
            while j < sorted.len() && sorted[j] == v {
                j += 1;
            }
            pts.push((v as f64, j as f64 / n));
            i = j;
        }
        pts
    }

    /// Plot one step-CDF per `n_init` to a PNG.
    pub fn plot_fig8_cdf(
        path: &str,
        series: &[(usize, Vec<(f64, f64)>)],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (900, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let x_max = series
            .iter()
            .flat_map(|(_, p)| p.iter().map(|&(x, _)| x))
            .fold(1.0_f64, f64::max);
        let mut chart = ChartBuilder::on(&root)
            .caption("Fig. 8 - Algorithm-2 iteration CDF", ("sans-serif", 28))
            .margin(12)
            .x_label_area_size(45)
            .y_label_area_size(55)
            .build_cartesian_2d(0f64..(x_max + 1.0), 0f64..1.02f64)?;
        chart
            .configure_mesh()
            .x_desc("iterations")
            .y_desc("empirical CDF")
            .draw()?;
        let palette = [RED, BLUE, GREEN];
        for (i, (n_init, pts)) in series.iter().enumerate() {
            let color = palette[i % palette.len()];
            chart
                .draw_series(LineSeries::new(pts.iter().cloned(), color.stroke_width(2)))?
                .label(format!("n_init = {n_init}"))
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], color));
        }
        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.85))
            .border_style(BLACK)
            .draw()?;
        root.present()?;
        Ok(())
    }

    /// Plot solve-time vs |T| on log-log axes to a PNG.
    pub fn plot_fig9_timing(
        path: &str,
        rows: &[Fig9Row],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (900, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let xs: Vec<f64> = rows.iter().map(|r| r.grid_len as f64).collect();
        let ys: Vec<f64> = rows.iter().map(|r| r.seconds.max(1e-6)).collect();
        let x_lo = xs.iter().cloned().fold(f64::INFINITY, f64::min).max(1.0);
        let x_hi = xs.iter().cloned().fold(0.0, f64::max).max(10.0);
        let y_lo = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_hi = ys.iter().cloned().fold(0.0, f64::max);
        let mut chart = ChartBuilder::on(&root)
            .caption("Fig. 9 - solve time vs |T|", ("sans-serif", 28))
            .margin(12)
            .x_label_area_size(45)
            .y_label_area_size(70)
            .build_cartesian_2d(
                (x_lo..x_hi * 1.5).log_scale(),
                (y_lo * 0.5..y_hi * 2.0).log_scale(),
            )?;
        chart
            .configure_mesh()
            .x_desc("|T| (grid size)")
            .y_desc("solve time [s]")
            .draw()?;
        chart.draw_series(LineSeries::new(
            xs.iter().cloned().zip(ys.iter().cloned()),
            BLUE.stroke_width(2),
        ))?;
        chart.draw_series(
            xs.iter()
                .cloned()
                .zip(ys.iter().cloned())
                .map(|(x, y)| Circle::new((x, y), 3, BLUE.filled())),
        )?;
        root.present()?;
        Ok(())
    }
```

- [ ] **Step 4: Wire PNG writes into `fig8` and `fig9`**

At the end of `fig8`, after the CSV write, add:

```rust
        let cdf: Vec<(usize, Vec<(f64, f64)>)> = N_INITS
            .iter()
            .map(|&n_init| {
                let counts: Vec<usize> = rows
                    .iter()
                    .filter(|r| r.n_init == n_init)
                    .map(|r| r.iterations)
                    .collect();
                (n_init, empirical_cdf(&counts))
            })
            .collect();
        match plot_fig8_cdf("target/fig8_cdf.png", &cdf) {
            Ok(()) => println!("  CDF plot             : target/fig8_cdf.png"),
            Err(e) => eprintln!("  PNG write failed     : {e}"),
        }
```

At the end of `fig9`, after the CSV write, add:

```rust
        match plot_fig9_timing("target/fig9_timing.png", &rows) {
            Ok(()) => println!("  timing plot          : target/fig9_timing.png"),
            Err(e) => eprintln!("  PNG write failed     : {e}"),
        }
```

- [ ] **Step 5: Run the CDF test and produce the PNGs**

Run: `cargo test --features validation --bin monte_carlo empirical_cdf -- --nocapture`
Expected: PASS.

Run: `cargo run --features validation --bin monte_carlo -- fig8`
Expected: also prints `CDF plot : target/fig8_cdf.png`; the file exists and is a nonempty PNG (`ls -l target/fig8_cdf.png`).

- [ ] **Step 6: Run the full gate**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green.

- [ ] **Step 7: Commit**

```bash
git add src/bin/monte_carlo.rs
git commit -m "feat(mc): optional plotters PNGs (iteration CDF + log-log timing)"
```

---

### Task 6: CI invariant integration test (`tests/monte_carlo.rs`)

**Files:**
- Create: `tests/monte_carlo.rs`

**Interfaces:**
- Consumes: the public API only (`solve`, `J2Roe`, `AbsoluteOrbit`, `Piecewise`, `Pseudostate`, `SolveParams`, `TimeGrid`) + `rand`/`rand_distr`. Samples independently of the bin (invariants hold for any seed).
- Produces: one `#[test] fn monte_carlo_invariants_hold()`.

- [ ] **Step 1: Write the integration test**

Create `tests/monte_carlo.rs`:

```rust
//! Phase 6 CI invariant test: Monte Carlo behaviour of the public `solve` API on the
//! worked-example problem. Asserts paper-INDEPENDENT invariants (NOT the paper's
//! 4.90/3.99/3.31 means — see spec §6 Phase 6 validation stance). Runs only under the
//! `validation` feature (for rand/rand_distr); CI runs `--all-features`.
#![cfg(feature = "validation")]

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_planner::{solve, Pseudostate, SolveParams, TimeGrid};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::TAU;

const A_C: f64 = 25_000e3;
const N_SAMPLES: usize = 64; // tunable for CI runtime (192 solves total)

fn chief() -> AbsoluteOrbit {
    AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    )
}

fn sample_ws(n: usize, seed: u64) -> Vec<Pseudostate> {
    let mut rng = StdRng::seed_from_u64(seed);
    let normal = Normal::new(0.0_f64, 1000.0).expect("σ > 0");
    (0..n)
        .map(|_| {
            let mut c = [0.0_f64; 6];
            for x in c.iter_mut() {
                *x = normal.sample(&mut rng);
            }
            Pseudostate::from_row_slice(&c) / A_C
        })
        .collect()
}

#[test]
fn monte_carlo_invariants_hold() {
    let dynamics = J2Roe::new(chief(), 0.0, 117_990.0);
    let cost = Piecewise::new(TAU / chief().mean_motion());
    let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0);
    let ws = sample_ws(N_SAMPLES, 0xC0FFEE);
    let n_inits = [2usize, 6, 10];

    let mut means = [0.0_f64; 3];
    let mut max_iters = 0usize;
    let mut max_res = 0.0_f64;
    let mut failures = 0usize;

    for (k, &n_init) in n_inits.iter().enumerate() {
        let params = SolveParams { n_init, ..SolveParams::default() };
        let (mut sum, mut count) = (0usize, 0usize);
        for &w in &ws {
            match solve(&dynamics, &cost, w, grid, &params) {
                Ok(sol) => {
                    sum += sol.iterations;
                    count += 1;
                    max_iters = max_iters.max(sol.iterations);
                    max_res = max_res.max(sol.residual);
                }
                Err(_) => failures += 1,
            }
        }
        assert!(count > 0, "n_init={n_init}: no successful solves");
        means[k] = sum as f64 / count as f64;
    }

    // Invariant 1: every solve succeeds (Phase 5b robustness).
    assert_eq!(failures, 0, "{failures} solve(s) failed; expected 0");
    // Invariant 2: converges within the paper's stated 8-iteration bound.
    assert!(max_iters <= 8, "max iterations {max_iters} exceeds the paper's 8-iter bound");
    // Invariant 3: residual under 0.01% (the min-fuel SOCP reconstructs w).
    assert!(max_res < 1e-4, "max residual {max_res:.3e} exceeds 0.01%");
    // Invariant 4: Fig. 8 shape — more init times ⇒ fewer iterations.
    assert!(
        means[0] > means[2],
        "mean iters n_init=2 ({:.3}) should exceed n_init=10 ({:.3})",
        means[0],
        means[2]
    );

    eprintln!(
        "observed mean iters: n_init=2 -> {:.2}, 6 -> {:.2}, 10 -> {:.2}  (paper 4.90/3.99/3.31)",
        means[0], means[1], means[2]
    );
}
```

- [ ] **Step 2: Run the test and observe the numbers**

Run: `cargo test --features validation --test monte_carlo -- --nocapture`
Expected: PASS, and the `eprintln!` reports the three observed mean iteration counts.

> **If an assertion fails (observe-then-lock, per the validation stance):**
> - `max_iters <= 8` fails → record the actual max; this is a real finding. First confirm it is not a degenerate-`w` artifact, then set the bound to the observed max and add a one-line spec note (Phase-5 band methodology). Do **not** silently widen without recording why.
> - `failures == 0` fails → investigate (likely an `n_init=2` coarse pick that doesn't span ⇒ unbounded dual SOCP ⇒ `SolverFailed`, the Phase-4 finding). Either document a small failure tolerance with the cause, or fix the seeding. Record the decision in the spec.
> - `means[0] > means[2]` fails → unexpected; investigate before relaxing.

- [ ] **Step 3: Run the full gate**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green, including `monte_carlo_invariants_hold` and the bin-internal unit tests.

- [ ] **Step 4: Commit**

```bash
git add tests/monte_carlo.rs
git commit -m "test(mc): seeded CI invariant test for Monte Carlo solve behaviour"
```

---

### Task 7: Finalize — spec update, full validation run, memory, handoff

**Files:**
- Modify: `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md` (mark Phase 6 done; record observed numbers)
- Modify: `memory/spec-validation-status.md` + `memory/MEMORY.md` (status update)

- [ ] **Step 1: Produce the full validation artifacts once**

Run: `cargo run --release --features validation --bin monte_carlo`
Expected: both sweeps run; `target/fig8_iterations.csv` (600 rows), `target/fig9_timing.csv` (6 rows), `target/fig8_cdf.png`, `target/fig9_timing.png` all written. Record from stdout: the three Fig. 8 mean iteration counts, `frac<=8`, max residual; and the Fig. 9 constant-then-linear shape (note the `10⁴`→`10⁶` rows).

- [ ] **Step 2: Update the spec's Phase 6 entry to ✅ Done**

In `docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`, change the `### Phase 6` heading to `✅ Done` and append a short results paragraph with the *observed* numbers from Step 1 (mean iters per `n_init` vs the paper's 4.90/3.99/3.31; whether all residuals were `< 0.01%`; the Fig. 9 shape; final test count). State explicitly whether the observed means landed near the paper's (and, if not, that the reframed stance applies — same as Phase 5). Update the top-of-file status line (`Phases 0–5b complete` → include Phase 6).

- [ ] **Step 3: Run the full gate one final time**

Run: `cargo fmt --all && cargo fmt --check && cargo clippy --all-features -- -D warnings && cargo build --all-features && cargo test --all-features`
Expected: all green; note the final total test count.

- [ ] **Step 4: Update memory**

Update `memory/spec-validation-status.md` (append a Phase 6 update: harness implemented, observed MC numbers, test count, branch/commit) and the matching one-line entry in `memory/MEMORY.md`. Convert any relative dates to absolute (today = 2026-06-19).

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md memory/
git commit -m "docs(spec): Phase 6 Monte Carlo harness done — observed MC results recorded"
```

- [ ] **Step 6: Hand off to branch completion**

Invoke the `superpowers:finishing-a-development-branch` skill to choose merge / PR / cleanup for `phase6-monte-carlo` (the established per-phase workflow: whole-branch review, then PR to `main`). Do **not** auto-merge.

---

## Self-Review

**Spec coverage** (each §6 Phase 6 requirement → task):
- Feature-gated bin + `#[cfg(not)]` stub → Task 1. ✓
- `rand`/`rand_distr` behind `validation`; portable `StdRng` + documented seed → Tasks 1, 2. ✓
- Sampling convention (`Normal(0, 1000 m)` ÷ `a_c`; never near-zero) → Task 2. ✓
- Fixed problem = Table III chief + eq. 49 cost → Task 1 builders, used everywhere. ✓
- Fig. 8: 200 paired `w` across `n_init ∈ {2,6,10}`, `n_coarse=20` fixed; record iters/residual/total_dv; count (don't panic) failures; report means vs 4.90/3.99/3.31 as reference; CSV; optional CDF PNG → Tasks 3, 5. ✓
- Fig. 9: sizes `10…10⁶`, `dt=(t_f−t_i)/(n−1)`, fixed Table III `w`, warmup + timed `Instant`, record actual `grid_len`; CSV; optional log-log PNG; documented 10⁶ cost; `fig8`/`fig9` arg → Tasks 4, 5. ✓
- CI invariant test: seeded, smaller `N`, asserts success / ≤8 iters / `<0.01%` residual / `mean(2)>mean(10)`; not the paper means; observe-then-lock → Task 6. ✓
- Determinism, serial (no rayon), exit criteria → Tasks 2/4/7. ✓

**Placeholder scan:** every code step shows complete code; every run step has an exact command + expected output. The only `#[allow(dead_code)]` is transient (Task 2 → removed Task 3, explicitly). The "observe-then-lock" note in Task 6 gives concrete default assertions plus a documented adjustment rule — not a placeholder. ✓

**Type consistency:** `Fig8Row`/`Fig8Stat`/`Fig9Row` field names match between definition (Tasks 3/4), the CSV writers, the summaries, and the plot wiring (Task 5). `run_fig8`/`summarize_fig8`/`write_fig8_csv`/`fig8`, `run_fig9`/`write_fig9_csv`/`fig9`, `empirical_cdf`/`plot_fig8_cdf`/`plot_fig9_timing`, `sample_pseudostates` are spelled identically at definition and call sites. All planner API names match `mdot.rs`/`types.rs` (`J2Roe::new`, `AbsoluteOrbit::new`, `mean_motion`, `Piecewise::new`, `solve`, `Solution.{iterations,residual,total_dv}`, `TimeGrid::uniform`/`.len`, `SolveParams{n_init,..}`, `Pseudostate::from_row_slice`). ✓
