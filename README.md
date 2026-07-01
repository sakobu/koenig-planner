# koenig-damico-planner

[![crates.io](https://img.shields.io/crates/v/koenig-damico-planner.svg)](https://crates.io/crates/koenig-damico-planner)
[![docs.rs](https://img.shields.io/docsrs/koenig-damico-planner)](https://docs.rs/koenig-damico-planner)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/koenig-damico-planner.svg)](#license)

> Minimum-Δv impulsive maneuver planning for spacecraft relative orbits — a
> finite-difference-verified Rust port of the Koenig & D'Amico (IEEE TAC 2020)
> fuel-optimal impulsive control algorithm.

```sh
cargo add koenig-damico-planner
```

A faithful Rust re-implementation of Koenig & D'Amico's fuel-optimal impulsive
control algorithm for linear systems with time-varying cost
(_"Fast Algorithm for Fuel-Optimal Impulsive Control of Linear Systems with
Time-Varying Cost,"_ IEEE Transactions on Automatic Control, 2020).

It computes a minimum-Δv impulsive maneuver plan that drives a deputy spacecraft
to a target set of quasi-nonsingular relative orbital elements (ROEs) under J2
secular dynamics, via the paper's three-step reachable-set method:

1. **Initialize** a candidate time grid (Algorithm 1),
2. **Refine** it by solving the dual reachability SOCP and adding/dropping contact times until convergence (Algorithm 2),
3. **Extract** the maneuvers by a direct gauge-aware minimum-fuel SOCP over the converged active set (Algorithm 3).

## Status & fidelity

- **Dynamics are finite-difference verified.** The J2 mean-ROE state-transition
  matrix `Φ(t,t_f)` and control-input matrix `B(t)` are independently checked by
  `tests/fd_stm.rs` and `tests/fd_b_matrix.rs` at two orbit regimes.
- **STM correction.** The paper's printed `Φ₂₁` δλ-drift term is a transcription
  typo (`−1.5 n Δt²`, dimensionally invalid); this crate uses the correct linear
  `−1.5 n Δt` (the printed Δt² form is dimensionally inconsistent for a rate term). See `src/dynamics/stm.rs`.
- **The worked-example figures are not bit-reproducible.** Under the corrected
  dynamics the paper's §VIII Table IV maneuvers do not reconstruct the Table III
  target, which is consistent with transcription errors in the published example. The crate validates
  the _math_ and _self-consistency_, not the printed numbers — see
  `tests/worked_example.rs`.

## Usage

```rust,ignore
use koenig_damico_planner::{solve, SolveParams, TimeGrid};
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::cost::Piecewise;
use koenig_damico_planner::Pseudostate;
use std::f64::consts::TAU;

let a_c = 25_000e3; // chief semimajor axis [m], the I/O scale factor
let chief = AbsoluteOrbit::new(
    a_c, 0.7, 40f64.to_radians(), 358f64.to_radians(), 0.0, 180f64.to_radians(),
);
let dynamics = J2Roe::new(chief, 0.0, 117_990.0)?;        // fallible: validates the chief
let grid = TimeGrid::uniform(0.0, 117_990.0, 30.0)?;       // fallible: validates dt>0, t_f>t_i
let cost = Piecewise::new(TAU / chief.mean_motion())?;     // fallible: validates period > 0
let w = Pseudostate::from_row_slice(&[50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]) / a_c;

let solution = solve(&dynamics, &cost, w, grid, &SolveParams::default())?;
println!("{} maneuvers, total dv = {:.4} mm/s",
    solution.maneuvers.len(), solution.total_dv * 1e3);
# Ok::<(), koenig_damico_planner::PlannerError>(())
```

A runnable version is `examples/mdot.rs`:

```sh
cargo run --example mdot
```

## Python bindings

The solver is also callable from Python (via [PyO3](https://pyo3.rs)), running the same native Rust
code locally — nothing is sent anywhere:

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install maturin
maturin develop -m crates/py/Cargo.toml      # build + install into the venv
```

```python
import koenig_planner as kp

chief = kp.Orbit(a=25_000e3, e=0.7, i=40.0, raan=358.0, argp=0.0, mean_anom=180.0)
#   a [m]; i, raan, argp, mean_anom in DEGREES.
sol = kp.solve(chief, t_i=0.0, t_f=117_990.0, dt=30.0,
               w_meters=[50, 5000, 100, 100, 0, 400], cost="piecewise")
print(sol.total_dv, "m/s in", len(sol.maneuvers), "maneuvers")
```

The package ships PEP 561 type stubs (`py.typed` + `.pyi`) for full editor/`mypy` support. See
[`crates/py/README.md`](crates/py/README.md) for details.

## HTTP service

The solver is also exposed as a self-hostable HTTP service (axum) — the same native code, running on
your own machine, nothing sent anywhere:

```bash
cargo run -p koenig-damico-planner-server   # listens on 0.0.0.0:8080 (override with KOENIG_PLANNER_ADDR)
```

```bash
curl -s localhost:8080/health               # {"status":"ok"}
curl -s -H 'content-type: application/json' \
     -d @crates/server/golden.json \
     localhost:8080/solve                    # POST a SolveRequest, get a SolveResponse
```

A `cargo-chef` → distroless **Dockerfile** is included for containerised self-hosting. See
[`crates/server/README.md`](crates/server/README.md) for the endpoints, the `{kind, message}` error
contract, and Docker usage.

## WASM browser demo

The solver also runs **entirely in the browser** via WebAssembly — the same native code compiled to
`wasm32`, solving client-side with nothing sent anywhere. The bindings are fully typed: generated
TypeScript types ([tsify](https://github.com/madonoharu/tsify)) give a `solve(req: SolveRequest): SolveOutcome`
that never throws — the error is returned as a value, so it stays typed in TypeScript — plus a
`solve_json(json): string` escape hatch.

```bash
wasm-pack build crates/wasm --target web     # → crates/wasm/pkg/ (wasm + generated .d.ts)
cd crates/wasm/www && npm install && npm run dev
```

The included demo (React + React-Three-Fiber, built with Vite) takes a chief orbit, target ROEs, and a
cost model and visualizes the plan client-side: two interactive 3D scenes — the chief orbit in ECI
and the deputy's relative orbit in RTN (chief at origin), both with per-maneuver Δv arrows (positioned
via the core's own Kepler solver) and a swept primer arrow — a playback scrubber, a Δv timeline, per-maneuver RTN
Δv components, and primer-vector panels (the primer magnitude vs time with the `|p| = 1` optimality
bound, plus the RTN primer components). See [`crates/wasm/README.md`](crates/wasm/README.md).

## Workspace layout

This repository is a Cargo workspace. The core solver is the root crate; the others are thin
frontends over a shared serde/JSON facade:

| Crate                                | Path                | Distribution                                                | Purpose                                                               |
| ------------------------------------ | ------------------- | ----------------------------------------------------------- | --------------------------------------------------------------------- |
| `koenig-damico-planner`              | `.` (root)          | [crates.io](https://crates.io/crates/koenig-damico-planner) | the core solver (this README)                                         |
| `koenig-damico-planner-api`          | `crates/api`        | internal (`publish = false`)                                | shared serde/JSON facade — the one `run()` / `run_json()` entry point |
| `koenig-damico-planner-py`           | `crates/py`         | planned — PyPI `koenig-planner` (import `koenig_planner`)   | Python bindings (above)                                               |
| `koenig-damico-planner-server`       | `crates/server`     | internal (`publish = false`)                                | self-hostable HTTP service (axum) — `POST /solve`, `GET /health`      |
| `koenig-damico-planner-wasm`         | `crates/wasm`       | npm, as `koenig-planner` (crates.io `publish = false`)      | WASM bindings + in-browser demo — `tsify`-typed `solve` / `solve_json` |
| `koenig-damico-planner-validation`   | `crates/validation` | internal (`publish = false`)                                | Monte-Carlo Fig. 7/8/9 reproduction harness                           |

## Validation harness (Fig. 7 / Fig. 8 / Fig. 9)

A Monte-Carlo harness in the workspace-internal `koenig-damico-planner-validation` crate
reproduces the paper's Fig. 7 (contact curve), Fig. 8 (refinement-iteration CDF under three
seedings), and Fig. 9 (solve time vs. grid size). Figure/CSV output is behind its `figures`
feature (which pulls `plotters`/`csv`):

    cargo run --release -p koenig-damico-planner-validation --features figures --bin monte_carlo

The seeded invariant test (`crates/validation/tests/invariants.rs`) asserts only
paper-independent invariants — the paper's reported means are shown as a _reference_, not a
pass/fail target.

**Fig. 8 — Algorithm-2 iteration distribution.** Empirical CDF of refinement
iterations over 200 random targets per seeding scheme (red n=2 window endpoints;
blue n=6 largest-contact times, i.e. Algorithm 1; green n=10 evenly spaced). Every
solve converges in ≤ 7 iterations — well within the paper's stated 8-iteration
bound — and mean iterations are 4.29 / 3.64 / 2.95 vs. the paper's 4.90 / 3.99 /
3.31.

![Fig. 8 — empirical CDF of Algorithm-2 refinement iterations for the three seeding schemes (n=2, n=6, n=10)](https://raw.githubusercontent.com/sakobu/koenig-planner/main/assets/fig8_cdf.png)

**Fig. 9 — solve time vs. grid size.** Wall-clock solve time for the Table III
target as the candidate-time grid grows from 10 to 10⁶ points (log–log, release
build). Solve time grows slowly at small |T| (setup-dominated), then scales
roughly linearly in |T| (~0.3 s at 10⁶ points).

![Fig. 9 — log–log plot of solve time versus candidate-time grid size, from 10 to 1e6 points](https://raw.githubusercontent.com/sakobu/koenig-planner/main/assets/fig9_timing.png)

## Build & test

```sh
cargo test --all-targets --all-features
```

## Stability & MSRV

This is a pre-1.0 (0.x) crate: breaking changes are batched into minor `0.x`
bumps and called out in [`CHANGELOG.md`](CHANGELOG.md); the serialized JSON wire
format is versioned with the crate. The public API exposes
[nalgebra](https://docs.rs/nalgebra) types, so a `nalgebra` major bump is a
breaking change of this crate. **MSRV: Rust 1.92**, enforced in CI and raised
only in a minor release.

## License

Licensed under either of [Apache-2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT) at
your option.
