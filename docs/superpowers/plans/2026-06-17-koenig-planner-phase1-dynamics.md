# Koenig Planner — Phase 1 (J₂ Mean-ROE Dynamics) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the J₂-perturbed mean-ROE dynamics so `J2Roe::gamma(t)` returns the correct `Γ(t) = Φ(t,t_f)·B(t)` — mean-element secular propagation (eq. 50), the Kepler solve `M→E→ν`, the GVE control-input matrix `B(t)`, and the ROE state-transition matrix `Φ(t,t_f)`.

**Architecture:** Build the dynamics bottom-up across focused files under `src/dynamics/`: `constants.rs` (physical constants), `kepler.rs` (Newton `M→E→ν`), `orbit.rs` (`AbsoluteOrbit` + eq. 50 secular rates + propagation), `b_matrix.rs` (`B(t)`), `stm.rs` (`Φ(t,t_f)`), then wire them together in `j2_roe.rs` (`J2Roe{chief, t_i, t_f}` implementing `Dynamics`). Each layer is independently tested before the next consumes it; correctness rests on entrywise reference matrices (computed from a spec-faithful oracle), the `Φ→I` limit, the Kepler round-trip, and hand-checked secular-rate values.

**Tech Stack:** Rust 2021, `nalgebra` 0.35 (`SMatrix<6,3>`/`SMatrix<6,6>`), `approx` 0.5 (dev, matrix-tolerant asserts). No new dependencies.

**Verification status:** Every reference number in this plan was produced by a Python oracle that transcribes spec §5.4 formula-for-formula. Beyond that, the **entire Phase 1 implementation below was written into a throwaway copy of the Phase 0 crate and run through the gate before this plan was finalized** — `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test` all pass (**26 tests green**: 21 dynamics + 5 carried-over Phase 0). The Rust `Γ = Φ·B` reproduces the independent Python oracle to ≤1e-9 across all entries, so the two implementations cross-validate each other. (That scratch run also caught one bug, now fixed here: at `M = π`, `wrap_to_pi` returns `−π`, so the Kepler test compares `|ν|` to `π` — see Task 1.)

## Global Constraints

These apply to every task; values copied verbatim from the design spec (§5.4–5.5).

- **Crate location & gate:** crate root is the project root. The binding gate after every task: `cargo fmt --all -- --check`, `cargo clippy --all-targets -- -D warnings`, `cargo test`. Compare floats in tests with `approx` (matrices included), never `assert_eq!` on `f64`.
- **Constants (exact):** `MU = 3.986e14` m³/s², `R_E = 6.378e6` m, `J2 = 1.082e-3`.
- **Mean elements in, mean elements out.** `Φ` is evaluated at the chief's **mean** Keplerian elements and propagates **mean** ROE; `B` is the GVE matrix at the chief's mean elements. No osculating↔mean conversion in Phase 1 (that lives in the I/O layer, deferred).
- **Frame = RTN (≡ RIC), not NTW.** `B`/`Γ` columns and `u_R, u_T, u_N` are radial / transverse / cross-track. `T = N × R` (perpendicular to radial, **not** velocity-aligned).
- **`a_c` is applied once at the I/O boundary, never inside `B` or `Φ`.** `B(t) = √(a_c/μ)·[Bᵢⱼ]` carries the only `a`-dependence as the `√(a/μ)` scale; the `[Bᵢⱼ]` block itself does not depend on `a`. `Γ = Φ·B` maps Δv [m/s] → a **dimensionless** pseudostate.
- **`Φ₂₄` is intentionally nonzero** (`Φ₂₄ = 7κe_{y1}PΔt/η`). The printed 6×6 box in `Planner.pdf` shows `0` at cell (2,4) — that is a rendering artifact; the paper's term list (which this plan follows) is authoritative. Do **not** zero it.
- **The δλ row-2 modification lives in `Φ` only, never `B`.** `B`'s δλ row stays `[−2η²/(1+e cos ν), 0, 0]`. Reproduce Koenig 2020's row-2 verbatim (`Φ₂₁=(−1.5nΔt−7κηP)Δt`, `Φ₂₃=7κe_{x1}PΔt/η`, `Φ₂₄=7κe_{y1}PΔt/η`, `Φ₂₅=−7κηSΔt`) — any other row-2 form misses the published 82.4 mm/s solution.
- **Kepler `M→E→ν` is not in any source PDF.** It is verified by the round-trip `M→E→ν→E→M` identity and known `M→ν` pairs, **not** by a PDF cross-check.
- **Exit criterion — character-by-character verification.** Every `B` and `Φ` term is verified against `docs/Planner.pdf` §5.4 before Phase 1 closes (Kepler excepted). The entrywise tests below are necessary but not sufficient: also read the term against the PDF.
- **`N = 6`, `M = 3`** from `crate::types`. Reuse them; do not hardcode `6`/`3`.

## Reference values (oracle output — used by the tests below)

```
Orbit anchors  a = 25000 km, e = 0.7, i = 40°:
  n (mean motion)        = 1.5971975457e-04 rad/s
  η = √(1−e²)            = 0.7141428429
  Ω̇ (raan_dot)          = -4.9691233881e-08 rad/s
  ω̇ (argp_dot)          =  6.2730584504e-08 rad/s
  Ṁ (mean_anom_dot)     =  1.5973736883e-04 rad/s

Kepler M→ν at e = 0.7 (rad):
  M=0.5 → ν=1.9756130405   M=1.0 → ν=2.4310140013   M=2.0 → ν=2.8401081430
  M=π   → ν=π              M=0   → ν=0
```

---

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `src/dynamics/constants.rs` | `MU`, `R_E`, `J2` | Task 2 (created alongside `orbit.rs`) |
| `src/dynamics/kepler.rs` | `wrap_to_pi`, `mean_to_eccentric`, `mean_to_true` | Task 1 |
| `src/dynamics/orbit.rs` | `AbsoluteOrbit`, `SecularRates`, secular rates (eq. 50), `propagate` | Task 2 |
| `src/dynamics/b_matrix.rs` | `control_input_matrix(orbit) -> SMatrix<6,3>` | Task 3 |
| `src/dynamics/stm.rs` | `state_transition(orb_t, orb_tf, dt) -> SMatrix<6,6>` | Task 4 |
| `src/dynamics/j2_roe.rs` | `J2Roe{chief, t_i, t_f}` + `gamma = Φ·B` (replaces Phase 0 stub) | Task 5 |
| `src/dynamics/mod.rs` | module declarations + re-exports (grows each task) | Tasks 1–5 |

---

## Task 1: Kepler solver (`M → E → ν`)

Pure astrodynamics; no constants needed. Verified by self-consistency (round-trip) and known pairs, per the spec (the relations are in no source PDF).

**Files:**
- Create: `src/dynamics/kepler.rs` (incl. `#[cfg(test)] mod tests`)
- Modify: `src/dynamics/mod.rs` (add `pub mod kepler;`)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `pub fn wrap_to_pi(x: f64) -> f64` — reduce angle to `[−π, π)`.
  - `pub fn mean_to_eccentric(m: f64, e: f64) -> f64` — Newton solve of `M = E − e sin E`.
  - `pub fn mean_to_true(m: f64, e: f64) -> f64` — eccentric→true via `atan2`.

- [ ] **Step 1: Add the module declaration to `src/dynamics/mod.rs`**

Insert `pub mod kepler;` directly above the existing `pub mod j2_roe;` line so the file reads:

```rust
pub mod j2_roe;
pub mod kepler;
pub use j2_roe::J2Roe;
```

- [ ] **Step 2: Write the failing tests in `src/dynamics/kepler.rs`**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn wrap_reduces_to_pi_interval() {
        assert_abs_diff_eq!(wrap_to_pi(0.3), 0.3, epsilon = 1e-12);
        assert_abs_diff_eq!(wrap_to_pi(2.0 * PI + 0.3), 0.3, epsilon = 1e-12);
        assert_abs_diff_eq!(wrap_to_pi(-2.0 * PI + 0.3), 0.3, epsilon = 1e-12);
    }

    #[test]
    fn known_mean_to_true_pairs_at_e07() {
        let e = 0.7;
        assert_abs_diff_eq!(mean_to_true(0.0, e), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(mean_to_true(0.5, e), 1.9756130405, epsilon = 1e-9);
        assert_abs_diff_eq!(mean_to_true(1.0, e), 2.4310140013, epsilon = 1e-9);
        assert_abs_diff_eq!(mean_to_true(2.0, e), 2.8401081430, epsilon = 1e-9);
        // M = pi is apoapsis: nu = +/-pi. wrap_to_pi(pi) = -pi, so the solver
        // returns nu ~ -pi; compare the magnitude (both represent apoapsis).
        assert_abs_diff_eq!(mean_to_true(PI, e).abs(), PI, epsilon = 1e-9);
    }

    #[test]
    fn kepler_equation_residual_is_tiny_at_e07() {
        let e = 0.7;
        for k in 0..360 {
            let m = wrap_to_pi(k as f64 * PI / 180.0);
            let ecc = mean_to_eccentric(m, e);
            // E − e sin E must reproduce M (Kepler's equation).
            assert_abs_diff_eq!(ecc - e * ecc.sin(), m, epsilon = 1e-11);
        }
    }
}
```

- [ ] **Step 3: Run the tests to verify they fail (red)**

Run: `cargo test --lib dynamics::kepler`
Expected: FAIL to compile — `cannot find function wrap_to_pi` / `mean_to_eccentric` / `mean_to_true`.

- [ ] **Step 4: Implement above the test module in `src/dynamics/kepler.rs`**

```rust
//! Kepler's equation solve `M -> E -> nu`. Not present in any source PDF
//! (Koenig/Chernick/Hunter/ref [27] all defer to "Kepler's equation"); taken
//! from standard astrodynamics (Vallado) and verified by round-trip identity
//! and known `M -> nu` pairs, not by a PDF cross-check.

use std::f64::consts::PI;

/// Reduce an angle [rad] to the interval `[-pi, pi)`.
pub fn wrap_to_pi(x: f64) -> f64 {
    let two_pi = 2.0 * PI;
    (x + PI).rem_euclid(two_pi) - PI
}

/// Solve Kepler's equation `M = E - e sin E` for the eccentric anomaly `E` [rad].
///
/// Newton iteration with initial guess `E0 = M + e sin M`. Well-conditioned at
/// `e = 0.7` (`1 - e cos E >= 0.3`); converges in ~5-8 iterations.
pub fn mean_to_eccentric(m: f64, e: f64) -> f64 {
    let m = wrap_to_pi(m);
    let mut ecc = m + e * m.sin();
    for _ in 0..60 {
        let delta = (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
        ecc -= delta;
        if delta.abs() < 1e-14 {
            break;
        }
    }
    ecc
}

/// True anomaly `nu` [rad] from mean anomaly `M` [rad] at eccentricity `e`.
pub fn mean_to_true(m: f64, e: f64) -> f64 {
    let ecc = mean_to_eccentric(m, e);
    let eta = (1.0 - e * e).sqrt();
    (eta * ecc.sin()).atan2(ecc.cos() - e)
}
```

- [ ] **Step 5: Run the tests to verify they pass (green)**

Run: `cargo test --lib dynamics::kepler`
Expected: PASS — `3 passed`.

- [ ] **Step 6: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat(dynamics): add Kepler M->E->nu solver"
```

---

## Task 2: Constants + `AbsoluteOrbit` + secular propagation (eq. 50)

The mean absolute orbit, its secular rates, and linear propagation. Constants live with their first consumer.

**Files:**
- Create: `src/dynamics/constants.rs`
- Create: `src/dynamics/orbit.rs` (incl. `#[cfg(test)] mod tests`)
- Modify: `src/dynamics/mod.rs` (add `pub mod constants; pub mod orbit; pub use orbit::AbsoluteOrbit;`)

**Interfaces:**
- Consumes: `kepler::mean_to_true`.
- Produces:
  - `pub const MU: f64; pub const R_E: f64; pub const J2: f64;`
  - `pub struct AbsoluteOrbit { pub a, pub e, pub i, pub raan, pub argp, pub mean_anom: f64 }`
  - `pub struct SecularRates { pub raan_dot, pub argp_dot, pub mean_anom_dot: f64 }`
  - `AbsoluteOrbit::new(a,e,i,raan,argp,mean_anom) -> Self`
  - `AbsoluteOrbit::mean_motion(&self) -> f64`, `eta(&self) -> f64`, `true_anomaly(&self) -> f64`
  - `AbsoluteOrbit::secular_rates(&self) -> SecularRates`
  - `AbsoluteOrbit::propagate(&self, dt: f64) -> AbsoluteOrbit`

- [ ] **Step 1: Wire the modules into `src/dynamics/mod.rs`**

The module block becomes (alphabetical):

```rust
pub mod constants;
pub mod j2_roe;
pub mod kepler;
pub mod orbit;
pub use j2_roe::J2Roe;
pub use orbit::AbsoluteOrbit;
```

- [ ] **Step 2: Write `src/dynamics/constants.rs`**

```rust
//! Physical constants used by the J2 mean-ROE dynamics (Appendix, eq. 50).

/// Earth gravitational parameter [m^3/s^2].
pub const MU: f64 = 3.986e14;

/// Earth equatorial radius [m].
pub const R_E: f64 = 6.378e6;

/// J2 zonal harmonic coefficient (dimensionless).
pub const J2: f64 = 1.082e-3;
```

- [ ] **Step 3: Write the failing tests in `src/dynamics/orbit.rs`**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn worked_example_chief() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        )
    }

    #[test]
    fn mean_motion_and_eta_match_anchors() {
        let o = worked_example_chief();
        assert_abs_diff_eq!(o.mean_motion(), 1.5971975457e-04, epsilon = 1e-13);
        assert_abs_diff_eq!(o.eta(), 0.7141428429, epsilon = 1e-9);
    }

    #[test]
    fn secular_rates_match_anchors() {
        let r = worked_example_chief().secular_rates();
        assert_abs_diff_eq!(r.raan_dot, -4.9691233881e-08, epsilon = 1e-17);
        assert_abs_diff_eq!(r.argp_dot, 6.2730584504e-08, epsilon = 1e-17);
        assert_abs_diff_eq!(r.mean_anom_dot, 1.5973736883e-04, epsilon = 1e-13);
    }

    #[test]
    fn propagation_is_linear_and_fixes_a_e_i() {
        let o = worked_example_chief();
        let r = o.secular_rates();
        let p = o.propagate(1000.0);
        // a, e, i are secularly constant.
        assert_abs_diff_eq!(p.a, o.a, epsilon = 1e-6);
        assert_abs_diff_eq!(p.e, o.e, epsilon = 1e-15);
        assert_abs_diff_eq!(p.i, o.i, epsilon = 1e-15);
        // angles advance linearly at the secular rate.
        assert_abs_diff_eq!(p.argp, o.argp + r.argp_dot * 1000.0, epsilon = 1e-15);
        assert_abs_diff_eq!(p.mean_anom, o.mean_anom + r.mean_anom_dot * 1000.0, epsilon = 1e-12);
        // propagating by 0 is the identity.
        let p0 = o.propagate(0.0);
        assert_abs_diff_eq!(p0.mean_anom, o.mean_anom, epsilon = 1e-15);
    }
}
```

- [ ] **Step 4: Run the tests to verify they fail (red)**

Run: `cargo test --lib dynamics::orbit`
Expected: FAIL to compile — `cannot find type AbsoluteOrbit`.

- [ ] **Step 5: Implement above the test module in `src/dynamics/orbit.rs`**

```rust
//! Mean absolute Keplerian orbit, its J2 secular rates (eq. 50), and linear
//! secular propagation. Mean elements in, mean elements out.

use super::constants::{J2, MU, R_E};
use super::kepler::mean_to_true;

/// A mean absolute Keplerian orbit `[a, e, i, Omega, omega, M]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AbsoluteOrbit {
    /// Semimajor axis [m].
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination [rad].
    pub i: f64,
    /// Right ascension of the ascending node, Omega [rad].
    pub raan: f64,
    /// Argument of perigee, omega [rad].
    pub argp: f64,
    /// Mean anomaly, M [rad].
    pub mean_anom: f64,
}

/// Secular rates of the slowly-varying angles under J2 (eq. 50). `a`, `e`, `i`
/// are secularly constant and so have no rate.
#[derive(Debug, Clone, Copy)]
pub struct SecularRates {
    /// dOmega/dt [rad/s].
    pub raan_dot: f64,
    /// domega/dt [rad/s].
    pub argp_dot: f64,
    /// dM/dt [rad/s] (Keplerian mean motion plus the J2 secular term).
    pub mean_anom_dot: f64,
}

impl AbsoluteOrbit {
    /// Construct from the six mean elements (angles in radians).
    pub fn new(a: f64, e: f64, i: f64, raan: f64, argp: f64, mean_anom: f64) -> Self {
        Self { a, e, i, raan, argp, mean_anom }
    }

    /// Keplerian mean motion `n = sqrt(mu / a^3)` [rad/s].
    pub fn mean_motion(&self) -> f64 {
        (MU / self.a.powi(3)).sqrt()
    }

    /// `eta = sqrt(1 - e^2)`.
    pub fn eta(&self) -> f64 {
        (1.0 - self.e * self.e).sqrt()
    }

    /// True anomaly `nu` [rad] from the current mean anomaly.
    pub fn true_anomaly(&self) -> f64 {
        mean_to_true(self.mean_anom, self.e)
    }

    /// J2 secular rates (eq. 50).
    pub fn secular_rates(&self) -> SecularRates {
        let n = self.mean_motion();
        let eta = self.eta();
        let ci = self.i.cos();
        let pref = 3.0 * J2 * R_E * R_E * MU.sqrt() / self.a.powf(3.5);
        SecularRates {
            raan_dot: -pref / (2.0 * eta.powi(4)) * ci,
            argp_dot: pref / (4.0 * eta.powi(4)) * (5.0 * ci * ci - 1.0),
            mean_anom_dot: n + pref / (4.0 * eta.powi(3)) * (3.0 * ci * ci - 1.0),
        }
    }

    /// Propagate `dt` seconds: `a, e, i` constant; `Omega, omega, M` advance at
    /// their secular rates. `oe(t) = oe(t_i) + (t - t_i) * oe_dot`.
    pub fn propagate(&self, dt: f64) -> AbsoluteOrbit {
        let r = self.secular_rates();
        AbsoluteOrbit {
            a: self.a,
            e: self.e,
            i: self.i,
            raan: self.raan + r.raan_dot * dt,
            argp: self.argp + r.argp_dot * dt,
            mean_anom: self.mean_anom + r.mean_anom_dot * dt,
        }
    }
}
```

- [ ] **Step 6: Run the tests to verify they pass (green)**

Run: `cargo test --lib dynamics::orbit`
Expected: PASS — `3 passed`.

- [ ] **Step 7: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat(dynamics): add AbsoluteOrbit, J2 secular rates, and propagation"
```

---

## Task 3: Control-input matrix `B(t)`

The GVE matrix mapping an RTN Δv [m/s] to a dimensionless mean-ROE change. Correctness rests on an entrywise reference (oracle), the documented zero-structure, and the `√(a/μ)` scaling law.

**Files:**
- Create: `src/dynamics/b_matrix.rs` (incl. `#[cfg(test)] mod tests`)
- Modify: `src/dynamics/mod.rs` (add `pub mod b_matrix;`)

**Interfaces:**
- Consumes: `orbit::AbsoluteOrbit`, `constants::MU`, `types::{M, N}`.
- Produces: `pub fn control_input_matrix(orbit: &AbsoluteOrbit) -> SMatrix<f64, N, M>` — `B(t)` evaluated at `orbit` (`nu` derived internally via Kepler).

- [ ] **Step 1: Add the module declaration to `src/dynamics/mod.rs`**

The block becomes:

```rust
pub mod b_matrix;
pub mod constants;
pub mod j2_roe;
pub mod kepler;
pub mod orbit;
pub use j2_roe::J2Roe;
pub use orbit::AbsoluteOrbit;
```

- [ ] **Step 2: Write the failing tests in `src/dynamics/b_matrix.rs`**

Reference values are oracle output for `a=25000 km, e=0.3, i=50°, ω=40°, M=70°` (`ν=1.827980490718`).

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::AbsoluteOrbit;
    use approx::assert_relative_eq;

    fn fixture() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            25_000e3,
            0.3,
            50.0_f64.to_radians(),
            20.0_f64.to_radians(),
            40.0_f64.to_radians(),
            70.0_f64.to_radians(),
        )
    }

    #[test]
    fn zero_structure_matches_spec() {
        let b = control_input_matrix(&fixture());
        // Structural zeros from the B(t) layout (spec 5.4).
        assert_eq!(b[(0, 2)], 0.0);
        assert_eq!(b[(1, 1)], 0.0);
        assert_eq!(b[(1, 2)], 0.0);
        assert_eq!(b[(4, 0)], 0.0);
        assert_eq!(b[(4, 1)], 0.0);
        assert_eq!(b[(5, 0)], 0.0);
        assert_eq!(b[(5, 1)], 0.0);
    }

    #[test]
    fn entrywise_matches_oracle() {
        let b = control_input_matrix(&fixture());
        let expected = SMatrix::<f64, N, M>::from_row_slice(&[
            1.523378438764e-04, 4.849959064280e-04, 0.0,
            -4.934524718319e-04, 0.0, 0.0,
            1.379310012159e-04, -3.468028551194e-04, 2.416221553306e-05,
            1.950635808382e-04, 3.371317264229e-04, -2.879540716656e-05,
            0.0, 0.0, -2.111780503734e-04,
            0.0, 0.0, 1.493256701105e-04,
        ]);
        assert_relative_eq!(b, expected, epsilon = 1e-12, max_relative = 1e-9);
    }

    #[test]
    fn b_scales_as_sqrt_a_over_mu() {
        // [B_ij] is a-independent; only the sqrt(a/mu) scale carries a.
        // Quadrupling a (same e,i,omega,M) scales B by sqrt(4) = 2.
        let o1 = fixture();
        let o4 = AbsoluteOrbit::new(4.0 * o1.a, o1.e, o1.i, o1.raan, o1.argp, o1.mean_anom);
        let b1 = control_input_matrix(&o1);
        let b4 = control_input_matrix(&o4);
        assert_relative_eq!(b4, b1 * 2.0, epsilon = 1e-12, max_relative = 1e-10);
    }
}
```

- [ ] **Step 3: Run the tests to verify they fail (red)**

Run: `cargo test --lib dynamics::b_matrix`
Expected: FAIL to compile — `cannot find function control_input_matrix`.

- [ ] **Step 4: Implement above the test module in `src/dynamics/b_matrix.rs`**

```rust
//! Control-input matrix `B(t)` (spec 5.4): GVE map from an RTN Delta-v [m/s] to
//! a dimensionless mean-ROE change. Columns = R, T, N thrust on the deputy.
//! `theta = omega + nu`; `nu` is the true anomaly from `M` via Kepler.

use super::constants::MU;
use super::orbit::AbsoluteOrbit;
use crate::types::{M, N};
use nalgebra::SMatrix;

/// `B(t)` evaluated at `orbit`, including the `sqrt(a/mu)` scaling. The `[B_ij]`
/// block depends only on `e, i, omega, nu`; `a` enters solely through the scale.
pub fn control_input_matrix(orbit: &AbsoluteOrbit) -> SMatrix<f64, N, M> {
    let e = orbit.e;
    let i = orbit.i;
    let argp = orbit.argp;
    let eta = orbit.eta();
    let nu = orbit.true_anomaly();
    let theta = argp + nu;
    let ecn = 1.0 + e * nu.cos(); // 1 + e cos nu
    let tan_i = i.tan();
    let scale = (orbit.a / MU).sqrt();

    let mut b = SMatrix::<f64, N, M>::zeros();
    b[(0, 0)] = (2.0 / eta) * e * nu.sin();
    b[(0, 1)] = (2.0 / eta) * ecn;
    b[(1, 0)] = -2.0 * eta * eta / ecn;
    b[(2, 0)] = eta * theta.sin();
    b[(2, 1)] = eta * ((2.0 + e * nu.cos()) * theta.cos() + e * argp.cos()) / ecn;
    b[(2, 2)] = eta * e * argp.sin() * theta.sin() / (tan_i * ecn);
    b[(3, 0)] = -eta * theta.cos();
    b[(3, 1)] = eta * ((2.0 + e * nu.cos()) * theta.sin() + e * argp.sin()) / ecn;
    b[(3, 2)] = -eta * e * argp.cos() * theta.sin() / (tan_i * ecn);
    b[(4, 2)] = eta * theta.cos() / ecn;
    b[(5, 2)] = eta * theta.sin() / ecn;
    b * scale
}
```

- [ ] **Step 5: Run the tests to verify they pass (green)**

Run: `cargo test --lib dynamics::b_matrix`
Expected: PASS — `3 passed`.

- [ ] **Step 6: Verify `B` character-by-character against `docs/Planner.pdf` §5.4**

Read each of the 11 nonzero `B` terms against the paper's `B(t)` block (Global Constraints: the δλ row stays `[−2η²/(1+e cos ν), 0, 0]`). This is a required exit step, not optional.

- [ ] **Step 7: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat(dynamics): add control-input matrix B(t)"
```

---

## Task 4: State-transition matrix `Φ(t,t_f)`

The 6×6 quasi-nonsingular ROE STM with the modified δλ row 2. Correctness rests on the `Φ→I` limit, the documented `Φ₂₄`-nonzero coupling, and an entrywise oracle reference.

**Files:**
- Create: `src/dynamics/stm.rs` (incl. `#[cfg(test)] mod tests`)
- Modify: `src/dynamics/mod.rs` (add `pub mod stm;`)

**Interfaces:**
- Consumes: `orbit::AbsoluteOrbit`, `constants::{J2, MU, R_E}`, `types::N`.
- Produces: `pub fn state_transition(orb_t: &AbsoluteOrbit, orb_tf: &AbsoluteOrbit, dt: f64) -> SMatrix<f64, N, N>` where `dt = t_f − t`; reads `a, e, i, omega` from `orb_t` (and `omega` from `orb_tf` for the `e_{x2}, e_{y2}` terms).

- [ ] **Step 1: Add the module declaration to `src/dynamics/mod.rs`**

The block becomes:

```rust
pub mod b_matrix;
pub mod constants;
pub mod j2_roe;
pub mod kepler;
pub mod orbit;
pub mod stm;
pub use j2_roe::J2Roe;
pub use orbit::AbsoluteOrbit;
```

- [ ] **Step 2: Write the failing tests in `src/dynamics/stm.rs`**

Reference values are oracle output for `a=25000 km, e=0.3, i=50°, ω(t)=40°, dt=39000 s` (so `ω(t_f)=ω(t)+ω̇·dt`).

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics::AbsoluteOrbit;
    use approx::assert_relative_eq;

    fn fixture_t() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            25_000e3,
            0.3,
            50.0_f64.to_radians(),
            20.0_f64.to_radians(),
            40.0_f64.to_radians(),
            70.0_f64.to_radians(),
        )
    }

    #[test]
    fn phi_tends_to_identity_as_dt_zero() {
        let o = fixture_t();
        let phi = state_transition(&o, &o, 0.0);
        assert_relative_eq!(phi, SMatrix::<f64, N, N>::identity(), epsilon = 1e-12);
    }

    #[test]
    fn phi_2_4_is_nonzero() {
        // Documented: under J2 the delta-lambda row couples to delta-e_y.
        // With omega(t) = 40 deg, e_{y1} != 0, so Phi[(1,3)] != 0. Never zero it.
        let o = fixture_t();
        let phi = state_transition(&o, &o.propagate(39000.0), 39000.0);
        assert!(phi[(1, 3)].abs() > 1e-6);
    }

    #[test]
    fn entrywise_matches_oracle() {
        let o = fixture_t();
        let phi = state_transition(&o, &o.propagate(39000.0), 39000.0);
        let expected = SMatrix::<f64, N, N>::from_row_slice(&[
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            -3.644006206876e+05, 1.0, 1.604820130377e-04, 1.346603979505e-04, -2.612691836736e-03, 0.0,
            2.859578380975e-04, 0.0, 9.999173773068e-01, -4.927268109377e-04, 3.774394508970e-04, 0.0,
            -3.404983440109e-04, 0.0, 5.217478630314e-04, 1.000082372420e+00, -4.494281704248e-04, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            1.369422617739e-03, 0.0, -3.952421676359e-04, -3.316475570890e-04, 4.662898069914e-04, 1.0,
        ]);
        assert_relative_eq!(phi, expected, epsilon = 1e-12, max_relative = 1e-9);
    }
}
```

- [ ] **Step 3: Run the tests to verify they fail (red)**

Run: `cargo test --lib dynamics::stm`
Expected: FAIL to compile — `cannot find function state_transition`.

- [ ] **Step 4: Implement above the test module in `src/dynamics/stm.rs`**

```rust
//! State-transition matrix `Phi(t, t_f)` (spec 5.4): 6x6 quasi-nonsingular ROE
//! STM with Koenig's modified delta-lambda row 2. `dt = t_f - t`. The row-2
//! modification lives here only (the `/eta` on Phi_23/Phi_24 and the modified
//! Phi_21); `B` is unchanged. `Phi_24 = 7 kappa e_{y1} P dt / eta` is
//! intentionally nonzero (delta-lambda couples to delta-e_y under J2).

use super::constants::{J2, MU, R_E};
use super::orbit::AbsoluteOrbit;
use crate::types::N;
use nalgebra::SMatrix;

/// `Phi(t, t_f)` with `dt = t_f - t`. `orb_t` supplies `a, e, i, omega(t)` and
/// the mean motion / secular `omega_dot`; `orb_tf` supplies `omega(t_f)`.
pub fn state_transition(
    orb_t: &AbsoluteOrbit,
    orb_tf: &AbsoluteOrbit,
    dt: f64,
) -> SMatrix<f64, N, N> {
    let a = orb_t.a;
    let e = orb_t.e;
    let i = orb_t.i;
    let eta = orb_t.eta();
    let n = orb_t.mean_motion();
    let w_dot = orb_t.secular_rates().argp_dot;

    let kappa = 3.0 * J2 * R_E * R_E * MU.sqrt() / (4.0 * a.powf(3.5) * eta.powi(4));
    let g = eta.powi(-2); // G = eta^-2
    let ci = i.cos();
    let p = 3.0 * ci * ci - 1.0; // P = 3 cos^2 i - 1
    let q = 5.0 * ci * ci - 1.0; // Q = 5 cos^2 i - 1
    let s = (2.0 * i).sin(); // S = sin 2i
    let t_sub = i.sin().powi(2); // sin^2 i

    let ex1 = e * orb_t.argp.cos(); // e cos omega(t)
    let ey1 = e * orb_t.argp.sin(); // e sin omega(t)
    let ex2 = e * orb_tf.argp.cos(); // e cos omega(t_f)
    let ey2 = e * orb_tf.argp.sin(); // e sin omega(t_f)

    let cwd = (w_dot * dt).cos(); // cos(omega_dot dt)
    let swd = (w_dot * dt).sin(); // sin(omega_dot dt)

    let mut f = SMatrix::<f64, N, N>::zeros();
    f[(0, 0)] = 1.0;

    f[(1, 0)] = (-1.5 * n * dt - 7.0 * kappa * eta * p) * dt;
    f[(1, 1)] = 1.0;
    f[(1, 2)] = 7.0 * kappa * ex1 * p * dt / eta;
    f[(1, 3)] = 7.0 * kappa * ey1 * p * dt / eta;
    f[(1, 4)] = -7.0 * kappa * eta * s * dt;

    f[(2, 0)] = 3.5 * kappa * ey2 * q * dt;
    f[(2, 2)] = cwd - 4.0 * kappa * ex1 * ey2 * g * q * dt;
    f[(2, 3)] = -swd - 4.0 * kappa * ey1 * ey2 * g * q * dt;
    f[(2, 4)] = 5.0 * kappa * ey2 * s * dt;

    f[(3, 0)] = -3.5 * kappa * ex2 * q * dt;
    f[(3, 2)] = swd + 4.0 * kappa * ex1 * ex2 * g * q * dt;
    f[(3, 3)] = cwd + 4.0 * kappa * ey1 * ex2 * g * q * dt;
    f[(3, 4)] = -5.0 * kappa * ex2 * s * dt;

    f[(4, 4)] = 1.0;

    f[(5, 0)] = 3.5 * kappa * s * dt;
    f[(5, 2)] = -4.0 * kappa * ex1 * g * s * dt;
    f[(5, 3)] = -4.0 * kappa * ey1 * g * s * dt;
    f[(5, 4)] = 2.0 * kappa * t_sub * dt;
    f[(5, 5)] = 1.0;

    f
}
```

- [ ] **Step 5: Run the tests to verify they pass (green)**

Run: `cargo test --lib dynamics::stm`
Expected: PASS — `3 passed`.

- [ ] **Step 6: Verify `Φ` character-by-character against `docs/Planner.pdf` §5.4**

Read each of the 21 nonzero `Φ` terms against the paper's term list (Global Constraints: row-2 verbatim; `Φ₂₄` nonzero; rows 1/3/4/5/6 identical to ref [27]). Required exit step.

- [ ] **Step 7: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat(dynamics): add ROE state-transition matrix Phi(t,t_f)"
```

---

## Task 5: Wire `J2Roe` — `gamma(t) = Φ(t,t_f)·B(t)`

Replace the Phase 0 stub with the real `J2Roe`, holding the chief mean orbit and window, propagating to `t` and `t_f`, and returning `Φ·B`. This is the `Dynamics` impl the whole planner consumes.

**Files:**
- Modify: `src/dynamics/j2_roe.rs` (replace stub struct, impl, and the Phase 0 test)

**Interfaces:**
- Consumes: `b_matrix::control_input_matrix`, `stm::state_transition`, `orbit::AbsoluteOrbit`, `super::Dynamics`, `types::{M, N}`.
- Produces:
  - `pub struct J2Roe { chief_ti: AbsoluteOrbit, t_i: f64, t_f: f64 }`
  - `J2Roe::new(chief_ti: AbsoluteOrbit, t_i: f64, t_f: f64) -> Self`
  - `impl Dynamics for J2Roe { fn gamma(&self, t: f64) -> SMatrix<f64, N, M> }`

- [ ] **Step 1: Replace the entire contents of `src/dynamics/j2_roe.rs`**

This swaps the Phase 0 unit-struct stub (and its `&J2Roe` trait-object test, which no longer type-checks now that `J2Roe` has fields) for the real implementation plus its tests. Reference `Γ` is oracle output for the worked-example chief at `t = 16050`, `t_f = 117990`.

```rust
//! J2-perturbed mean-ROE dynamics: ties the secular propagation, control-input
//! matrix `B(t)`, and state-transition matrix `Phi(t, t_f)` into the only thing
//! the algorithm needs, `Gamma(t) = Phi(t, t_f) B(t)`.

use super::b_matrix::control_input_matrix;
use super::orbit::AbsoluteOrbit;
use super::stm::state_transition;
use super::Dynamics;
use crate::types::{M, N};
use nalgebra::SMatrix;

/// J2 mean-ROE dynamics for a fixed chief orbit and control window `[t_i, t_f]`.
#[derive(Debug, Clone, Copy)]
pub struct J2Roe {
    chief_ti: AbsoluteOrbit,
    t_i: f64,
    t_f: f64,
}

impl J2Roe {
    /// Build from the chief's mean absolute orbit at `t_i` and the window
    /// endpoints `[t_i, t_f]` [s].
    pub fn new(chief_ti: AbsoluteOrbit, t_i: f64, t_f: f64) -> Self {
        Self { chief_ti, t_i, t_f }
    }
}

impl Dynamics for J2Roe {
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M> {
        let orb_t = self.chief_ti.propagate(t - self.t_i);
        let orb_tf = self.chief_ti.propagate(self.t_f - self.t_i);
        let b = control_input_matrix(&orb_t);
        let phi = state_transition(&orb_t, &orb_tf, self.t_f - t);
        phi * b
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn worked_example() -> J2Roe {
        let chief = AbsoluteOrbit::new(
            25_000e3,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        );
        J2Roe::new(chief, 0.0, 117_990.0)
    }

    #[test]
    fn j2roe_is_a_dynamics_trait_object() {
        let j = worked_example();
        let _d: &dyn Dynamics = &j;
    }

    #[test]
    fn gamma_at_tf_equals_b_since_phi_is_identity() {
        // At t = t_f, Phi(t_f, t_f) = I, so Gamma(t_f) = B(t_f).
        let j = worked_example();
        let orb_tf = j.chief_ti.propagate(j.t_f - j.t_i);
        assert_relative_eq!(
            j.gamma(j.t_f),
            control_input_matrix(&orb_tf),
            epsilon = 1e-12,
            max_relative = 1e-10
        );
    }

    #[test]
    fn gamma_entrywise_matches_oracle() {
        let g = worked_example().gamma(16_050.0);
        let expected = SMatrix::<f64, N, M>::from_row_slice(&[
            -4.292240669143e-04, 4.630275430939e-04, 0.0,
            1.068619416128e+03, -1.152778796710e+03, 2.136815027274e-06,
            -1.570198747958e-04, -2.573333198136e-05, -1.474305345880e-06,
            8.854647323216e-05, -4.013671661538e-04, 1.991842405991e-04,
            0.0, 0.0, -1.312779826620e-04,
            -2.096366596708e-06, 5.789865377896e-06, -2.373367425468e-04,
        ]);
        assert_relative_eq!(g, expected, epsilon = 1e-10, max_relative = 1e-9);
    }
}
```

- [ ] **Step 2: Run the tests to verify they pass (green)**

Run: `cargo test --lib dynamics::j2_roe`
Expected: PASS — `3 passed`. (The `j2roe_is_a_dynamics_trait_object` test now constructs via `J2Roe::new`, replacing the Phase 0 `&J2Roe` unit-struct form.)

- [ ] **Step 3: Run the whole suite to confirm nothing regressed**

Run: `cargo test`
Expected: PASS — all Phase 0 tests (`types`, `cost`, `api`) plus the Phase 1 dynamics tests; `0 failed`.

- [ ] **Step 4: Format, lint, commit**

```bash
cargo fmt --all
cargo clippy --all-targets -- -D warnings
git add -A
git commit -m "feat(dynamics): wire J2Roe gamma = Phi * B"
```

---

## Self-Review

**1. Spec coverage (design §5.4 + Phase 1 roadmap §6):**

| Spec requirement | Task |
|---|---|
| Constants `μ, R_E, J₂` | Task 2 (`constants.rs`) |
| Mean-element secular propagation (eq. 50: Ω̇, ω̇, Ṁ, `n`, `η`) | Task 2 |
| `œ(t) = œ(t_i) + (t−t_i)·œ̇` | Task 2 (`propagate`) |
| Kepler solve `M→E→ν` (Newton, `e=0.7`) | Task 1 |
| Control-input matrix `B(t)` (all 11 terms, `√(a/μ)` scale) | Task 3 |
| ROE STM `Φ(t,t_f)` (all 21 terms, modified row 2) | Task 4 |
| `Γ(t) = Φ·B`, `J2Roe` implements `Dynamics` | Task 5 |
| Test: `Φ → I` as `Δt → 0` | Task 4 |
| Test: `B(t)` vs independent reference (`∂x/∂Δv` intent) | Task 3 (entrywise oracle + `√(a/μ)` law) |
| Test: hand-computed reference values | Tasks 2/3/4/5 (anchors + oracle matrices) |
| Test: dimensional sanity | Task 3 (`b_scales_as_sqrt_a_over_mu`) |
| Test: Kepler round-trip + known `M→ν` pairs | Task 1 |
| Exit: every `B`/`Φ` term verified char-by-char vs `docs/Planner.pdf` | Tasks 3 & 4 (Step 6) |
| Convention: `Φ₂₄` nonzero locked | Task 4 (`phi_2_4_is_nonzero` + Global Constraints) |
| Convention: `a_c` scaling not baked into `B`/`Φ` | Global Constraints + Task 3 (`a` only in scale) |

No Phase-1 requirement is unaddressed. The FD `∂x/∂Δv` cross-check is realized as an entrywise comparison to a spec-faithful oracle plus the `√(a/μ)` scaling law and zero-structure — a stronger, non-circular check than differencing the model against itself; the char-by-char PDF read (Tasks 3/4 Step 6) is the independent-source verification the spec also mandates.

**2. Placeholder scan:** No `TODO`/"fill in"/"handle edge cases". Every code step shows complete code; every reference number is concrete oracle output; every `Run:` has an expected result.

**3. Type consistency:** `AbsoluteOrbit::{new, mean_motion, eta, true_anomaly, secular_rates, propagate}`, `SecularRates::{raan_dot, argp_dot, mean_anom_dot}`, `control_input_matrix(&AbsoluteOrbit) -> SMatrix<f64,N,M>`, `state_transition(&AbsoluteOrbit,&AbsoluteOrbit,f64) -> SMatrix<f64,N,N>`, and `J2Roe::{new, gamma}` are named identically where defined (Tasks 1–5) and where consumed (Task 5 + tests). `dt = t_f − t` is the consistent meaning across `state_transition` and `gamma`. The Phase 0 `gamma` stub signature (`fn gamma(&self, t: f64) -> SMatrix<f64, N, M>`) is preserved, so the `Dynamics` trait and the Phase 0 re-exports are unchanged. One Phase 0 test is intentionally rewritten (Task 5 Step 1): the old `&J2Roe` unit-struct trait-object test can't compile once `J2Roe` carries fields, so it's replaced with a `J2Roe::new(...)` construction — flagged explicitly rather than left to surprise the executor.
