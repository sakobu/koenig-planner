# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — 2026-06-19

Initial release.

### Added
- Faithful Rust port of the Koenig–D'Amico fuel-optimal impulsive control
  algorithm (IEEE TAC 2020): the three-step reachable-set method — candidate
  time-grid initialization (Algorithm 1), dual-reachability SOCP refinement
  (Algorithm 2), and direct gauge-aware minimum-fuel SOCP maneuver extraction
  (Algorithm 3).
- J2 secular mean-ROE dynamics: state-transition matrix `Φ(t,t_f)` and
  control-input matrix `B(t)`, with a fallible Kepler→B→γ pipeline.
- Cost models: `Norm2`, `FaceMax`, and the perigee-windowed `Piecewise` cost.
- Public API: `solve`, `solve_from_initial_times`, and the convex-encoding
  building blocks (`min_fuel_socp`, `refine_socp`).
- Monte-Carlo validation harness (Fig. 8 / Fig. 9) behind the `validation`
  feature; a seeded, paper-independent invariant test runs in CI without it.
- Dual MIT / Apache-2.0 licensing.

### Fidelity notes
- **STM correction.** The paper's printed `Φ₂₁` δλ-drift term (`−1.5 n Δt²`) is
  a dimensionally-invalid transcription typo; this crate uses the correct linear
  `−1.5 n Δt`. The correction was confirmed by the paper's first author.
- **Finite-difference-verified dynamics.** `Φ(t,t_f)` and `B(t)` are checked by
  independent FD oracles (`tests/fd_stm.rs`, `tests/fd_b_matrix.rs`) at two orbit
  regimes.
- **Worked example is not bit-reproducible.** Under the corrected dynamics the
  paper's §VIII example does not reconstruct its own target, consistent with
  transcription errors in the published numbers; the crate validates the math
  and self-consistency rather than the printed figures.

[0.1.0]: https://github.com/sakobu/koenig-planner/releases/tag/v0.1.0
