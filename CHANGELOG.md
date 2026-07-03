# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `AbsoluteOrbit::time_to_perigee()` — the duration `[s]` from the orbit's epoch to
  its next perigee (`M ≡ 0`), referenced to the Keplerian mean motion. It is the
  single source of truth for placing the default eq.-49 piecewise perigee windows,
  reused by both the solver (`api::run`) and the WASM presentation geometry.
  Additive and non-breaking.

### Changed
- **Breaking (behavioral): planner output changes for `t_i ≠ 0` piecewise-default
  requests.** The default piecewise perigee epoch
  (`CostSpec::Piecewise { t_perigee0: None }`) now anchors the eq.-49 FaceMax
  windows at the chief's actual perigee in absolute grid time,
  `t_i + time_to_perigee()`. The previous default computed `(-M₀/n) mod period`,
  omitting `t_i`, so for a non-zero `t_i` the cheap-fuel perigee window was shifted
  by `t_i (mod period)` and the solver returned a valid but sub-optimal plan.
  Affects every frontend (all route through `api::run`). No public API changed
  (`cargo semver-checks` clean); this is a numerical-output change only, and it is
  **byte-identical for `t_i = 0`**, so the paper's worked example and every golden
  are unchanged. \[KD20\] eq. 49.
- **Breaking (behavioral): WASM presentation geometry changes for `t_i ≠ 0`.**
  `crates/wasm` `chief_geometry` now evaluates the chief/deputy at the `t_i` epoch
  by propagating each absolute grid time `t` as the duration `t - t_i`; previously
  it passed the absolute time directly to `AbsoluteOrbit::propagate` (which advances
  the epoch by a *duration*), so for a non-zero `t_i` every burn marker, playback
  track, and primer arrow was over-propagated by `t_i`. The `ChiefGeometry` values
  (npm `koenig-planner`) therefore change for `t_i ≠ 0` — presentation-only: the
  plan/maneuvers and the numerical core are unchanged.
- **Breaking (behavioral): WASM `solve` outcome and error kind.** A valid solve
  whose target ROE reconstructs a non-elliptic deputy (`e ≥ 1`) now returns `Ok`
  with the deputy-derived geometry fields (`deputy_track_rtn`, `maneuver_rtn`)
  empty, where it previously returned `Err { kind: "solver" }` and discarded the
  plan. Separately, a genuine (now-unreachable) presentation-geometry fault reports
  `kind: "internal"` instead of the mislabeled `"solver"`. Callers switching on the
  outcome `status` or `error.kind` for these cases will observe the new values; the
  `SolveOutcome` / `ApiErrorKind` types are unchanged.

## [0.4.0] — 2026-07-02

### Added
- WASM demo: the RTN relative-motion scene now draws per-maneuver burn markers
  with Δv (thrust) arrows and a swept primer arrow, matching the ECI scene (it
  previously showed only the relative orbit and the moving deputy glyph). The WASM
  `ChiefGeometry` gains two presentation fields — `maneuver_rtn` (a new
  `ManeuverRtnDto` carrying the deputy's relative burn position and the native-RTN
  Δv) and `primer_rtn` (the RTN primer history) — an additive, non-breaking change
  to the generated TypeScript types (npm `koenig-planner`). Presentation-only,
  reusing the core's Kepler solver; the numerical core is untouched.
- `InvalidInputKind::NonFiniteChiefAngle { name, value }` classifies a chief
  orbit whose right-ascension, argument of perigee, or mean anomaly is non-finite.
  Additive and non-breaking (`InvalidInputKind` is `#[non_exhaustive]`); it maps
  to `ErrorClass::InvalidInput` like every other input error.
- `InvalidInputKind::Tolerance { eps_cost, eps_remove }` classifies a refinement
  tolerance that is non-finite or `<= 0`. Additive and non-breaking
  (`InvalidInputKind` is `#[non_exhaustive]`); it maps to `ErrorClass::InvalidInput`
  like every other input error.

### Changed
- `dynamics::stm::state_transition` now takes the chief at `t` and `dt` only —
  `state_transition(orb_t, dt)` instead of `state_transition(orb_t, orb_tf, dt)`.
  It derives `ω(t_f) = ω(t) + ω_dot·dt` internally from the eq. 50 secular drift,
  so `e_x2`/`e_y2` and the `cos`/`sin(ω_dot·dt)` rotation now share one `ω_dot·dt`
  and a single mean orbit fully determines `Φ(t, t_f)`. The old three-argument form
  let a caller pass an `orb_tf` that was not `orb_t.propagate(dt)`: only
  `orb_tf.argp` was ever read, so a different `a`/`e`/`i` was silently ignored,
  yielding a matrix that is the STM of no trajectory. **Breaking** (public function
  signature) — the redundant argument is removed rather than runtime-checked, so the
  inconsistent state is now unrepresentable and `Φ` stays infallible. The per-entry
  STM and finite-difference oracles are byte-identical, and the worked-example
  summary is unchanged (residual 1.135e-14, total_dv 81.3542 mm/s); the live-path
  `Φ` can differ only at the f64-ULP level — a mathematically identical rerounding
  of `ω(t_f)` (the old path anchored it at `t_i` rather than deriving it over `dt`).

### Removed
- WASM `ChiefGeometry` no longer includes `relative_trajectory_rtn` (the deputy's
  relative orbit sampled over one chief period). It was recomputed on every
  `solve()` but never consumed by the demo, which draws the relative motion from
  `deputy_track_rtn` (sampled on the full-mission playback grid) instead. This is
  a **breaking** removal from the generated TypeScript types (npm `koenig-planner`);
  no runtime behavior changes beyond a smaller payload.

### Fixed
- `TimeGrid` no longer manufactures a candidate time past `t_f`. The grid-point
  count used `round((t_f - t_i)/dt) + 1`, so on a window whose length is not a
  whole multiple of `dt` the final sample `t_i + (len-1)·dt` could land up to
  `dt/2` **after** `t_f` — an inadmissible candidate (the paper restricts the
  control domain to `T ⊆ [t_i, t_f]`, eq. 5) that `J2Roe::gamma` then evaluated
  with a backward-extrapolated state-transition matrix (`t_f - t < 0`), so a
  maneuver could be scheduled after the epoch at which the target must already be
  achieved while the reported residual stayed near zero. The count is now
  floor-based (with a relative tolerance that preserves the exact endpoint on
  commensurate windows), so every candidate stays within `[t_i, t_f]`.
  **Behavior change:** on a non-commensurate window the grid now has one fewer
  point and its last sample falls short of `t_f` by less than `dt` (`t_f` is a
  grid point only when the window length is a whole multiple of `dt`); the
  `TimeGrid` docs are corrected to state this contract. The paper's worked example
  (`[0, 117990]` @ 30 s) and the Hunter cross-check (`[0, 39000]` @ 10 s) are
  commensurate and unaffected (solver output byte-identical; worked-example
  residual 1.135e-14).
- `J2Roe::new` now rejects a chief whose inclination, right-ascension, argument
  of perigee, or mean anomaly is non-finite (NaN or ±∞), returning
  `PlannerError::InvalidInput` (a caller-fixable "bad request"). Previously the
  inclination's `sin(i).abs() < 1e-9` singularity guard let a non-finite `i`
  through (`NaN < 1e-9` is false) and `Ω`/`ω`/`M` had no finiteness check at all,
  so a malformed chief either surfaced later as a mis-classified `KeplerDivergence`
  (`ErrorClass::Unsolvable`) once the NaN reached the Kepler solve, or — for `Ω`,
  which `B(t)` and the STM never read — produced a wrongly-successful solve.
  Rustdoc that described a validating `J2Roe` as making `KeplerDivergence`
  unreachable is corrected to the accurate guarantee (a malformed chief is
  rejected up front). No numerical change on valid inputs (worked-example residual
  1.135e-14).
- `solve` / `solve_from_initial_times` now validate `SolveParams.eps_cost` and
  `eps_remove`, rejecting any non-finite or non-positive tolerance as
  `PlannerError::InvalidInput(InvalidInputKind::Tolerance { .. })`. The paper's
  Algorithm 2 requires both strictly positive; previously they reached refinement
  unchecked, where a `NaN` tolerance burned all 50 refinement iterations and
  surfaced as a mis-classified `NotConverged` (`ErrorClass::Unsolvable`) with a
  `NaN` target (a negative tolerance did likewise) — a caller-fixable input now
  reported as one. No effect on valid runs (defaults `0.01` / `0.01`); no numerical
  change (worked-example residual 1.135e-14).
- The Kepler solver (`mean_to_eccentric`) is now globally convergent for every
  valid eccentricity `e ∈ [0, 1)`. The unglobalized Newton iteration diverged near
  periapsis for near-parabolic chiefs (`e >= ~0.995`; ~0.3 % of mean anomalies at
  `e = 0.999`), where the derivative `1 - e cos E` shrinks to ~1e-3 and a Newton
  step overshoots ~100×; a single such grid time failed the entire solve with
  `KeplerDivergence` (`ErrorClass::Unsolvable`), even though the chief was
  physically valid (e.g. a high-apogee `e = 0.995` orbit with perigee above
  Earth). The solve now falls back to bisection on the bracket `[-π, π]`, where
  `F(E) = E - e sin E - M` is strictly increasing (`F' = 1 - e cos E > 0`) with a
  single bracketed root, so convergence is guaranteed. The Newton fast path is
  unchanged, so every previously-converging input returns a bit-identical `E`
  (checked across 12M evaluations up to `e = 0.9999`); `KeplerDivergence` is
  retained as a now-unreachable regression backstop. Rustdoc that described the
  divergence as "not reachable for valid `e`" is corrected to the guarantee the
  solver now actually provides. No numerical change on valid inputs
  (worked-example residual 1.135e-14).

### Documentation
- WASM `ChiefGeometry` doc comments (which surface as JSDoc on the npm
  `koenig-planner` types) no longer call `perigee_window` / `perigee_arc_eci` the
  "FaceMax band". Per the paper (eq. 49), that band is the piecewise cost's
  perigee attitude-constraint window (T1), where the gauge switches to FaceMax
  (Norm2 elsewhere) — it is not produced by the standalone `FaceMax` cost. Doc
  wording only; no type or behavior change. The demo's ECI-scene comment was
  corrected to match.
- `J2Roe::new` documents the near-equatorial conditioning cliff in a new
  `# Conditioning` rustdoc section: the `sin(i)` guard rejects only the exact
  `1/tan(i)` singularity, but for a chief within `~1e-4` rad of equatorial the
  cross-track `B(t)` entries grow to `~1e5` (vs `~1e-4` typical), so the SOCP data
  is ill-conditioned and the solver may return degraded accuracy or a loud
  `SolverFailed` — never a silently wrong plan. Documents an existing numerical
  limitation; no behavior change.
- `solver::extract_qp` documents in a new `# Accuracy` section that it returns
  the raw Algorithm 3 QP optimum with no post-solve verification — the clarabel
  status is accepted at reduced accuracy (`AlmostSolved`) and the magnitudes are
  not residual- or budget-checked — matching the paper's gateless Algorithm 3 and
  unlike the gated live `solve` extraction. A nonzero weighted residual is the
  correct optimum when the budget or a nonnegativity bound is active, so callers
  needing a verified, exactly-reachable plan should use `solve`. Documents the
  existing contract of an off-`solve`-path primitive; no behavior change.
- `AbsoluteOrbit::propagate` documents in a new `# Angle range` section that
  `Ω, ω, M` are returned unbounded (faithful to eq. 50's linear form) and that
  this is safe because every consumer uses these angles only through `sin`/`cos`,
  with the Kepler solve re-wrapping `M` internally; callers needing a bounded
  angle are pointed at `wrap_to_pi`. No behavior change.

## [0.3.0] — 2026-06-27

> **Migrating from 0.2.0.** The error/gauge changes are breaking for direct Rust
> API consumers only; the `w_metres` → `w_meters` rename additionally breaks the
> JSON wire format, the Python keyword, and the WASM/TS request type.
>
> - `PlannerError` is now `#[non_exhaustive]` → add a trailing `_` arm to any
>   `match`, or classify with the new `PlannerError::class() -> ErrorClass`.
> - `SublevelSet` is now sealed → only the built-in `Norm2`/`FaceMax` gauges
>   implement it (this was never a documented extension point).
> - The request field `w_metres` is renamed to `w_meters` (US spelling) → rename
>   the key in any JSON request body, the Python `solve`/`solve_json` `w_meters=`
>   keyword, and the WASM/TS `SolveRequest` field. The field is required, so a
>   stale `w_metres` key fails loudly with a missing-field error.

### Added

- The crate now re-exports `nalgebra` (`koenig_damico_planner::nalgebra`) and
  documents that its public API exposes `nalgebra` types, so a `nalgebra` major
  bump is a breaking change of this crate (downstream can use the version-matched
  re-export).

### Changed

- **BREAKING:** `PlannerError` is now `#[non_exhaustive]`, so future error
  categories are non-breaking for direct Rust consumers. New public
  `PlannerError::class() -> ErrorClass` classifies an error into a coarse,
  transport-agnostic category (`InvalidInput` / `Unsolvable`); the api frontend
  maps it to the HTTP error kind, so a future core variant is still classified
  inside the core crate at compile time. `ErrorClass` is re-exported at the crate
  root. The serialized wire JSON is unchanged.
- **BREAKING:** `SublevelSet` is now a sealed trait — only the built-in
  `Norm2` / `FaceMax` gauges implement it — so future methods can be added to it
  without a breaking change. The seal is reversible (it can be opened in a minor
  release if a downstream gauge is requested). `CostModel` and `Dynamics` remain
  open extension points.
- **BREAKING (wire / Python / WASM):** the request field `w_metres` is renamed to
  `w_meters`, standardizing the codebase on US-English spelling. This changes the
  JSON request key, the Python `solve` / `solve_json` `w_meters=` keyword, and the
  WASM/TS `SolveRequest` field. The field is required (and the HTTP `SolveRequest`
  is `deny_unknown_fields`), so a stale `w_metres` key is rejected loudly rather
  than silently ignored. The Rust core's public API is unchanged.

### Fixed
- Stale docs refreshed to match the shipped 0.2.0 frontends (docs-only — no API,
  behavior, or wire change): the root README's WASM-demo section now describes the
  React + React-Three-Fiber 3D console (two interactive 3D scenes plus a playback
  scrubber) instead of the superseded flat SVG demo; the
  `koenig-damico-planner-api` README and manifest no longer call the HTTP and WASM
  frontends "planned" and correctly attribute `run_json` to the Python/WASM
  `solve_json` escape hatches (the typed `solve` and the HTTP server call `run`);
  the public `Maneuver.t` rustdoc documents it as an absolute grid time
  (`t_i + k·dt`), not "measured from `t_i`"; and the `solver` module doc lists all
  three solver wrappers and attributes Algorithm 3 to `min_fuel_socp` (not the
  legacy `extract_qp`).

### Documentation
- Workspace-wide documentation accuracy and completeness pass (docs-only — no
  API, behavior, or wire change). Corrected: the root README's
  `TimeGrid::uniform` precondition comment (`t_f > t_i`, not `t_f >= t_i`); the
  `primer_history` signature in this changelog (it returns
  `Result<PrimerHistory, PlannerError>`); the api README's
  `deny_unknown_fields` claim (the nested `cost` enum is exempt); and the
  validation crate's stale "multi-second" Fig. 9 timing note (the release-build
  10⁶-point solve is ~0.3 s, matching the README). Filled completeness gaps: a
  crate-level rustdoc quick-start, a `serde`-feature note, and concrete-type
  pointers on the landing page; a "Limits" note and the server-only `internal`
  error-kind clarification in the api README; the `SolveResponse` body, graceful
  shutdown, and full env-var list in the server README; the primer-vector /
  `iterations` / `residual` outputs, optional tuning args, and the
  `primer_history.py` example pointer in the Python README; the React +
  React-Three-Fiber demo description, `version()`, and `solve_json`'s
  throw-on-failure behavior in the wasm README; and doc comments on the
  validation result structs. The Python type stubs now declare the read-only
  PyO3 fields as read-only properties so type checkers flag illegal writes.

## [0.2.0] — 2026-06-24

> **Migrating from 0.1.0.** This is the first release with breaking changes for
> direct Rust API consumers. The serialized JSON wire format is unchanged.
>
> - The `validation` feature is removed → depend on the new
>   `koenig-damico-planner-validation` crate for the Monte-Carlo harness.
> - `Piecewise::new` / `Piecewise::with_perigee_epoch` now return `Result` →
>   handle or `unwrap` the result.
> - `koenig-damico-planner-api`'s `ApiError.kind` is a typed `ApiErrorKind` (was
>   `&'static str`) → match the enum.
> - `PlannerError::InvalidInput` now wraps `InvalidInputKind`, and the
>   `min_fuel_socp` / `refine_socp` / `MinFuelSolution` / `RefineSolution` /
>   `ConicRows` / `Dual` / `FuelGenerator` items moved off the crate root to the
>   `solver::` / `types::` paths → update `match` arms and `use` paths.

### Added
- Primer-vector history on every solve (the paper's Fig. 7): the new public
  `primer_history(dynamics, cost, grid, lambda) -> Result<PrimerHistory, PlannerError>` reconstructs the
  primer `p(t) = Γᵀ(t)·λ` and its dual-gauge magnitude `g_{U(1,t)}(p(t))` at each
  grid time from the converged dual. The HTTP/Python/WASM `SolveResponse` now
  carries three parallel, grid-aligned arrays — `primer_times`, `primer_magnitude`
  (dimensionless, `≤ 1 + eps_cost`, `≈ 1` at maneuver times), and `primer_rtn` (the
  RTN primer vector itself — not the executed thrust direction, which under the
  polytopic/piecewise gauge is the support image `s(Γᵀλ)`). The WASM demo plots
  both (magnitude-vs-time with the `|p| = 1` bound
  and per-burn markers, plus the RTN components), and a `crates/py` matplotlib
  example does the same. Touch-1-away-from-a-burn reveals plan flexibility. The
  core solve path and `Solution` are unchanged, so the Monte-Carlo harness is
  unaffected.
- Tight snapshot regressions for the worked-example solutions (Koenig Table III and the
  Hunter L2 cross-check): total Δv, residual ceiling, maneuver count, and per-maneuver
  times/magnitudes are now pinned alongside the existing paper-bound bands, so silent
  science drift fails the test instead of passing.
- HTTP server now catches handler/middleware panics via `CatchPanicLayer` and returns the uniform `{"kind":"internal"}` 500 (panic payload logged server-side, never sent to the client). Wire-enum tags are pinned by tests.
- New workspace-internal crate `koenig-damico-planner-validation` (`crates/validation`) holding
  the Monte-Carlo sampler, the Fig. 7/8/9 reproduction harness, and the seeded invariant test.
  Figure/CSV generation is behind its `figures` feature.

### Changed
- `ApiError.kind` is now a typed `ApiErrorKind` enum (was `&'static str`), matched exhaustively by every frontend. The serialized wire JSON is unchanged; this is a breaking change for direct Rust consumers of `koenig-damico-planner-api`.
- **BREAKING:** the `validation` feature is removed from `koenig-damico-planner`; the Monte-Carlo
  harness and its `rand`/`rand_distr`/`plotters`/`csv` dependencies move to the new
  `koenig-damico-planner-validation` crate. The published core now depends only on
  `nalgebra`, `clarabel`, `thiserror` (and optional `serde`).
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
- **BREAKING:** `PlannerError::InvalidInput` now wraps a typed `InvalidInputKind`
  enum (was an opaque `String`). The new public `InvalidInputKind` classifies the
  cause (grid, eccentricity, period, budget, …) and carries the offending
  value(s); it is re-exported at the crate root. Rust consumers that matched on or
  constructed `InvalidInput(String)` must update to the enum. The serialized wire
  JSON is unchanged.
- **BREAKING:** the crate-root re-export surface is trimmed. `min_fuel_socp`,
  `refine_socp`, `MinFuelSolution`, and `RefineSolution` are now reached via
  `koenig_damico_planner::solver::…`, and `ConicRows`, `Dual`, and `FuelGenerator`
  via `koenig_damico_planner::types::…`, instead of the crate root. Update `use`
  paths accordingly; the items themselves are unchanged.

### Fixed
- The Monte-Carlo Fig. 9 timing sweep now surfaces solver failures: `run_fig9` returns a
  failure count (matching `run_fig8`), and the driver warns and skips the timing plot when any
  solve fails, instead of silently plotting a NaN-bearing series.
- Stale docs: the README no longer claims the seeded invariant test "runs without the feature
  flag", and the api golden test comment no longer mislabels the dual lower bound (≈0.0808) as
  the total Δv (≈0.0814).
- The `piecewise` cost's default perigee-window epoch is now derived from the
  chief's mean anomaly `M₀` — the first perigee passage at or after `t = 0`,
  `(-M₀ / n) mod period` — instead of assuming the chief is at apogee at
  `t = 0`. This places the eq. 49 FaceMax windows on the correct orbital arc
  for any chief; it reduces exactly to `period / 2` for the worked example
  (`M₀ = 180°`), leaving that result unchanged. Adapter-only (`api` and `wasm`,
  inherited by the Python and HTTP frontends); the numerical core is untouched.
- `TimeGrid::uniform` now rejects a zero-length window (`t_f == t_i`), agreeing
  with `validate_inputs` and `J2Roe::new`, which already require `t_f > t_i`.
  No working solve path changes — the solver already rejected a single-point
  grid one step later; this only tightens the standalone constructor so a grid
  it accepts is never rejected downstream.
- `docs.rs` now renders the derived `serde` `Serialize`/`Deserialize` impls on the
  public wire types (`Solution`, `Maneuver`, `SolveParams`, `TimeGrid`,
  `PlannerError`, `AbsoluteOrbit`). The crate has no default feature, so docs.rs
  built with `serde` off and omitted these impls from the published docs; a
  `[package.metadata.docs.rs]` entry now enables the `serde` feature for the docs
  build only. Docs-render metadata only: no code, API, or feature change, and not
  semver-relevant; it takes effect on the next published version.
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
- The `run_json` library entrypoint — and the Python/WASM frontends built on it —
  now rejects request bodies larger than `MAX_REQUEST_BYTES` (1 MiB, exposed as a
  public constant) with `bad_request`, before any JSON parse or allocation. This
  caps the previously-unbounded body size on those entrypoints; the HTTP server
  already enforced a 64 KiB body limit.
- Request DTOs now use `#[serde(deny_unknown_fields)]` (`OrbitDto`, `SolveParamsDto`,
  `SolveRequest`), so unknown/typo'd fields are rejected as `bad_request` instead of
  silently ignored. This closes the wire-format "no `deny_unknown_fields`" hardening
  item and bounds the unknown-field skip path on the `run_json`/py/wasm entrypoints
  (defense-in-depth, alongside the `MAX_REQUEST_BYTES` cap above).
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

### Documentation
- The public fallible functions (`solve`, `solve_from_initial_times`, and the
  convex-encoding building blocks `extract_qp` / `min_fuel_socp` / `refine_socp`)
  now carry `# Errors` rustdoc listing the `PlannerError` variants each can return
  and when; a `missing_errors_doc` / `missing_panics_doc` lint enforces the
  convention. `PlannerError::InvalidInput` is documented as a caller-fixable
  "bad request — correct the inputs" signal whose wrapped `InvalidInputKind`
  classifies the cause.
- The request wire contract is documented more completely: `n_coarse` / `n_init`
  are marked inert when `initial_times` is supplied (that path bypasses
  Algorithm 1) on the api and wasm DTOs and the Python `solve` docstring/stub,
  and the api crate README gains a "Wire stability" note tying the JSON
  request/response shape to crate semver and naming the stable cost-model tags,
  error kinds, and field names.

## [0.1.0] — 2026-06-19

Initial release.

### Added
- Faithful Rust port of the Koenig-D'Amico fuel-optimal impulsive control
  algorithm (IEEE TAC 2020): the three-step reachable-set method — candidate
  time-grid initialization (Algorithm 1), dual-reachability SOCP refinement
  (Algorithm 2), and direct gauge-aware minimum-fuel SOCP maneuver extraction
  (Algorithm 3).
- J2 secular mean-ROE dynamics: state-transition matrix `Φ(t,t_f)` and
  control-input matrix `B(t)`, with a fallible Kepler→B→γ pipeline.
- Cost models: `Norm2`, `FaceMax`, and the perigee-windowed `Piecewise` cost.
- Public API: `solve`, `solve_from_initial_times`, and the convex-encoding
  building blocks (`extract_qp`, `min_fuel_socp`, `refine_socp`).
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

[Unreleased]: https://github.com/sakobu/koenig-planner/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/sakobu/koenig-planner/releases/tag/v0.3.0
[0.2.0]: https://github.com/sakobu/koenig-planner/releases/tag/v0.2.0
[0.1.0]: https://github.com/sakobu/koenig-planner/releases/tag/v0.1.0
