# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- HTTP server now catches handler/middleware panics via `CatchPanicLayer` and returns the uniform `{"kind":"internal"}` 500 (panic payload logged server-side, never sent to the client). Wire-enum tags are pinned by tests.

### Changed
- `ApiError.kind` is now a typed `ApiErrorKind` enum (was `&'static str`), matched exhaustively by every frontend. The serialized wire JSON is unchanged; this is a breaking change for direct Rust consumers of `koenig-damico-planner-api`.
- **BREAKING:** `Piecewise::new` and `Piecewise::with_perigee_epoch` now return
  `Result<Self, PlannerError>` instead of `Self`. They validate that `period` is
  finite and `> 0` and that the perigee epoch is finite (`InvalidInput`
  otherwise), matching the fallible `TimeGrid::uniform` / `J2Roe::new`
  constructors. This prevents a non-finite or non-positive period from silently
  corrupting the eq.-49 window selector: a zero/NaN period makes `in_perigee_window`
  `false` for every time (collapsing the cost to pure `Norm2`), and a negative
  period makes it `true` for every time (collapsing it to `FaceMax` everywhere).
  Callers must now handle or `unwrap` the result; the `api` and `wasm` adapters
  surface a bad period as a `bad_request`.

### Fixed
- `docs.rs` now renders the derived `serde` `Serialize`/`Deserialize` impls on the
  public wire types (`Solution`, `Maneuver`, `SolveParams`, `TimeGrid`,
  `PlannerError`, `AbsoluteOrbit`). The crate has no default feature, so docs.rs
  built with `serde` off and omitted these impls from the published docs; a
  `[package.metadata.docs.rs]` entry now enables the `serde` feature for the docs
  build only — targeted (not `all-features`) so the tooling-only `validation`
  native deps (`plotters`/fontconfig) are not pulled in. Docs-render metadata
  only: no code, API, or feature change, and not semver-relevant; it takes effect
  on the next published version.
- `J2Roe::new` now rejects a chief whose semimajor axis `a` is not finite and
  positive, returning `PlannerError::InvalidInput` — completing the
  bounded-ellipse precondition alongside the existing `e ∈ [0,1)` check
  (\[KD20\] eq. 50 needs `n = √(μ/a³)` and an `a^{7/2}` denominator, real and
  finite only for `a > 0`). Previously a non-positive or non-finite `a` passed
  the constructor and only surfaced downstream as a NaN-poisoned `Γ`, which the
  frontends mis-reported as a `solver` failure instead of the caller-fixable
  `bad_request` it is. The constructor signature is unchanged (already
  fallible), so callers using `?`/`unwrap` are unaffected.
- The two genuine self-consistency checks that were `debug_assert!` — and so were
  compiled out of the release binary — are now always on. Algorithm 3's
  primal/dual cross-check (the extracted min-fuel objective must agree with the
  refinement dual budget `c*` to within tolerance, by strong duality —
  \[KD20\] Theorems 1–3) now returns `PlannerError::SolverFailed` on a mismatch
  instead of vanishing in release, and the `Piecewise` period precondition is
  enforced (see Changed).
- `Solution.total_dv` now reports the minimized fuel-cost objective — the paper's
  "delta-v cost" `c*` (eq. 4): `Σ‖Δvⱼ‖₂` under the `Norm2` cost and the polytope
  gauge `Σθ` under `FaceMax` — instead of `Σ‖Δvⱼ‖₂` of the recovered net Δv. The
  two agree for `Norm2` and for single-vertex `FaceMax` burns; they diverge (the
  gauge is larger, by up to √3 per burn) only when a perigee-window maneuver
  combines ≥2 tetrahedron vertices, where the old value under-reported the
  optimized budget. The value is now measured on the full, pre-prune solution,
  consistent with `residual`. Affects the reported number in every frontend
  (core/api/server/wasm/py); maneuvers, residual, and dual are unchanged.
- Completed the Apache-2.0 license appendix copyright; the Python package version
  is now sourced solely from the crate manifest (`pyproject.toml` uses
  `dynamic = ["version"]`), so it can no longer drift from `crates/py/Cargo.toml`.

### Security
- CI now scans dependencies for security advisories, license compatibility, and
  source provenance with `cargo-deny`, on every push/PR and on a weekly schedule
  (the schedule re-checks the committed `Cargo.lock` for newly-published
  advisories against unchanged dependencies). The demo's npm tree is audited
  (non-blocking).
- All GitHub Actions and Docker base images are now pinned to immutable commit
  SHAs / digests and kept current by Dependabot, removing the floating-tag and
  rolling-branch (`@master`/`@stable`) supply-chain exposure.
- The server container image now runs as a non-root user
  (`gcr.io/distroless/cc-debian12:nonroot`, UID 65532).

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
