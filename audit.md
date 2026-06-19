# Audit: `koenig-planner` Rust crate

This report audits the `koenig-planner` Rust crate (~4670 LoC) — a self-described faithful re-implementation of Koenig & D'Amico's 2020 IEEE TAC fuel-optimal impulsive control algorithm — against its design spec (`docs/superpowers/specs/2026-06-17-koenig-planner-rust-design.md`, §5 math reference). The method was a 10-dimension multi-agent read of the source in full, with high-severity findings independently re-verified against the source.

## Executive summary

This is a high-quality research crate. The core science is sound: the J2 mean-ROE dynamics (STM Φ, control matrix B, Kepler solve, secular rates, Γ composition) are encoded term-by-term faithful to the spec and are independently guarded by two genuine finite-difference oracles that reconstruct Φ and B by mathematically distinct routes. The Phase-5 Φ₂₁ δλ-drift fix (the paper's dimensionally-invalid `−1.5nΔt²` typo → linear `−1.5nΔt`) is real, correctly applied (`src/dynamics/stm.rs:45`), and the kind of bug the FD guards now catch cleanly. The cost models form a correct gauge/support-function dual pair; the three-algorithm pipeline faithfully restructures the paper's do/while; and the convex-solver wrappers handle status non-silently with no `unwrap` on solver output. Cross-cutting code quality is unusually disciplined: zero `panic!`/`unsafe`/`todo!` in non-test code, a single well-structured `thiserror` enum, centralized input validation, and clippy `-D warnings` enforced in CI.

The audit surfaced **no correctness defects**. The findings are documentation, packaging, and defensive-programming gaps. The two Medium items are both hygiene: a missing README and missing license-text files for a crate whose entire value proposition is fidelity-with-caveats. The most consequential correctness-adjacent observation is documentation drift: spec §5 — explicitly headed "what the code must encode" — still describes the superseded Algorithm-3 QP rather than the Phase-5b min-fuel SOCP the code actually runs. Several Low/Info items recur across dimensions: unguarded domain singularities (Kepler `e≥1`, `B(t)` at `i=0/π`), a dead `KeplerDivergence` error variant, and the retained-but-superseded `extract_qp`/`SolveParams.q` surface.

### Severity scorecard (post-verification)

| Severity | Count |
|----------|-------|
| Critical | 0 |
| High     | 0 |
| Medium   | 2 |
| Low      | 18 |
| Info     | 19 |
| **Total** | **39** |

High-severity findings requiring adversarial verification: 0 (no Critical/High findings were raised). Verification counts: confirmed 0, partial 0, refuted 0.

## Build, lint & test gate

**Verdict: fully green.** Run locally on 2026-06-19 with the stable `cargo` toolchain, mirroring the CI pipeline in `.github/workflows/ci.yml`.

| Stage | Command | Result |
|-------|---------|--------|
| Format | `cargo fmt --all -- --check` | ✅ clean (exit 0) |
| Lint | `cargo clippy --all-targets --all-features -- -D warnings` | ✅ zero warnings (exit 0) |
| Tests | `cargo test --all-targets --all-features` | ✅ **115 / 115 passed** — 0 failed, 0 ignored |

Test breakdown: 89 lib unit + 5 `monte_carlo` bin + 9 `algorithm` + 5 `api` + 1 `fd_b_matrix` + 1 `fd_stm` + 1 `monte_carlo` invariant + 1 `solver` + 3 `worked_example` (the `mdot` example compiles; 0 test fns). The two independent finite-difference guards (`fd_stm`, `fd_b_matrix`), the worked-example regression (`worked_example`), and the seeded Monte-Carlo invariant test all pass; the Monte-Carlo invariant test is the slowest at ~12 s.

> Caveat (tracked as a Low finding below): `cargo doc --no-deps --all-features` emits **29 `unresolved link` warnings** — unit annotations such as `[s]`/`[m]`/`[rad]` parsed as intra-doc links. These do not affect the build/lint/test gate, which is green.

## Strengths

De-duplicated across dimensions, the genuinely strong aspects are:

- **Faithful, FD-verified dynamics.** Every nonzero STM term (`src/dynamics/stm.rs:43-67`) and all 13 `B(t)` terms (`src/dynamics/b_matrix.rs:24-34`) match spec §5.4 character-by-character, including the subtle points (η-modified `Φ₂₃`/`Φ₂₄`, intentionally-nonzero `Φ₂₄`, mixed initial/final ω subscripts, `√(a/μ)` scaling). The Phase-5 `Φ₂₁` linear-drift fix is correct and consistent with the spec.
- **Genuine independent oracles.** `tests/fd_stm.rs` and `tests/fd_b_matrix.rs` reconstruct Φ and B by mathematically distinct routes (eq.51 ROE-map Jacobians + secular propagation; RTN-perturbed COE round-trips) with tolerances tight enough to separate a real coefficient/sign bug from the FD noise floor — exactly the guard that would have caught the Phase-1 typo. Both deliberately exercise the Hunter chief (ω≠0) to activate couplings that vanish at the worked-example chief.
- **Correct convex geometry.** FaceMax's `contact` is numerically confirmed to be exactly the support function of `U=conv{0, V_vertex}`, and the `Polytope` fuel generator is the gauge of the same `U` — a verified dual pair. The code correctly routes around the paper's eq.48 `V_face` typo, using it only as a transcription cross-check.
- **Faithful, progress-guaranteed algorithm port.** Algorithm 2's convergence-before-drop/add ordering (`refine.rs:103-115`) makes `T^opt` exactly the active set that produced `λ_opt`; `violated_local_maxima` is plateau/endpoint-aware and proven to always re-add the global max, so the loop always progresses. Thresholds match §5.1 exactly.
- **Clean architecture & layering.** The module tree maps one-to-one onto spec §4; the algorithm layer never touches clarabel directly; `Dynamics::gamma` exposes only Γ(t); error types are structured and machine-inspectable.
- **Disciplined cross-cutting quality.** The core library (`src/` excluding the `bin/` harness) has no `panic!`/`unsafe`/`todo!`/`unwrap` on production paths — the only production fallibility is a single infallible documented `expect` (`solver/mod.rs:21`, "default clarabel settings are always valid"). The Monte-Carlo harness binary adds one runtime `assert!` plus one immediately-guarded `.unwrap()` (`bin/monte_carlo.rs:458–465`, on a vector the assert proves non-empty) and a documented `expect` (`σ > 0`). Beyond that: `Result`-based boundaries, centralized `validate_inputs` using `!is_finite()` guards, and numerical hazards handled by formulation (total `atan2`, guarded divide-by-norm). Both `panic!` sites in the tree are test-only.
- **Strong, honest test suite.** 115 tests; the "reported-not-asserted" validation stance is scientifically defensible and well-documented (`paper_table_iv_does_not_reconstruct` even encodes why the paper's figures don't reproduce). Determinism, solver status mapping, and the degenerate collinear-support case are all tested.
- **Faithful Monte-Carlo harness.** Correctly distinguishes the three Fig.8 seedings (not an `n_init` knob), uses portable deterministic `StdRng`, gates heavyweight deps behind `validation`, and asserts only paper-independent invariants.

## Findings by dimension

### Architecture & module design

The crate realizes the spec's "Approach A" faithfully with clean separation of concerns. The wrinkles are residue from the Phase-5b/6 evolution: a superseded extractor and an unused params field remain in the public surface, and the crate root re-exports some solver/cost encoding internals. None are correctness defects.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | Superseded `extract_qp` still `pub` and re-exported | `src/lib.rs:16`, `src/solver/mod.rs:4,8` | — |
| Low | `SolveParams.q` is a public, self-documented-dead field | `src/types.rs:81-84,93` | — |
| Info | Crate-root re-exports expose solver/cost implementation detail | `src/lib.rs:16-20` | — |
| Info | `SublevelSet` gained `fuel_generator` beyond spec §4.1 sketch | `src/cost/mod.rs:8-18` vs spec §4.1 | — |

- **`extract_qp` re-export (Low):** The fixed-direction QP is dead on the solve path (Algorithm 3 calls `min_fuel_socp`) but billed equally at the crate root with no supersession signal. Drop it from the root re-export (keep `solver::extract_qp` for its pinning test) or add a rustdoc note that `min_fuel_socp` is the production extractor.
- **`SolveParams.q` (Low):** A required-shaped public field documented as "Reserved; unused by the Phase-5b min-fuel extractor." Remove it now the QP is off the solve path, or mark it `#[deprecated]`-style so it doesn't read as a live tuning knob.

### Dynamics — J2 mean-ROE model

A careful, faithful encoding of spec §5.4 with no math defects found. The `Φ₂₁` fix and `Φ₂₄` coupling are correctly implemented; κ and the secular-rate prefactor are mutually consistent to machine precision; the FD tests are genuine guards. Findings are defensive-programming gaps around domain singularities, all out of the paper's bounded-elliptic scope.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | Kepler solver silently returns wrong-but-finite value for `e ≥ 1` | `src/dynamics/kepler.rs:18-29` | — |
| Low | `B(t)` divides by `tan(i)` with no guard (undefined at `i=0`, `180°`) | `src/dynamics/b_matrix.rs:20,29,32` | — |
| Info | Kepler Newton loop has no non-convergence signal | `src/dynamics/kepler.rs:21-28` | — |
| Info | Physical constants are low-precision (3-4 sig figs) | `src/dynamics/constants.rs:4-10` | — |

- **Kepler `e ≥ 1` (Low):** `mean_to_eccentric` uses the elliptic equation with no precondition; a hyperbolic input "converges" to a meaningless finite value with no NaN/panic/error, flowing silently into `B(t)`. Add `debug_assert!(e >= 0.0 && e < 1.0)` (and/or in `AbsoluteOrbit::new`), or document the elliptic-only precondition.
- **`B(t)` `tan(i)` singularity (Low):** `B₃₃`/`B₄₃` divide by `tan_i`; at equatorial inclinations these become `±inf`/NaN. This faithfully matches the spec's printed `1/tan i` form, so it is a guard/doc gap, not a transcription error. Document the `i ∈ (0, π)` precondition or `debug_assert` inclination is bounded away from 0 and π.

### Cost models (Table II, eq. 47–49)

A faithful, mathematically consistent, well-tested encoding. The gauge/support-function duality is numerically verified; the eq.48 `V_face` typo is correctly diagnosed and routed around. The substantive items are an API footgun and a micro-perf nit; neither affects correctness.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | `Piecewise` silently assumes chief-at-apogee (`M0=180°`) and `t_i=0` | `src/cost/piecewise.rs:7-39` | — |
| Info | `vertex_columns()` recomputed (8 sqrt + 4 allocs) per call | `src/cost/facemax.rs:18-27,30-65` | — |

- **`Piecewise` implicit phase (Low):** Window centers are hard-wired to `(k+0.5)*period` in absolute `t`, which lands on perigee only when the chief is at apogee at `t=0` and the time origin is `t_i=0`. All current tests satisfy both, but a caller with `M0≠180°` or nonzero `t_i` gets windows at the wrong true anomaly with no error. Add an explicit phase/epoch parameter (e.g. `with_phase(period, t_perigee0)`) or document and validate the precondition; at minimum state the `t_i=0` assumption in rustdoc.

### Convex solver wrappers (clarabel SOCP/QP)

Cone construction is correct and verified against the clarabel source; the min-fuel solver faithfully encodes the Phase-5b program; status handling is non-silent with no `unwrap` on output; input guards are thoughtful. Findings are acceptable, documented tradeoffs around `AlmostSolved` and a debug-only consistency check.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | `refine_socp` accepts `AlmostSolved` without re-checking dual feasibility | `src/solver/refine_socp.rs:78-82`, `src/solver/mod.rs:27-34` | — |
| Low | Min-fuel vs dual-budget self-consistency check is only a `debug_assert` | `src/algorithm/extract.rs:44-48` | — |
| Info | Wrappers use clarabel default iter cap/tolerances (no headroom for large degenerate sets) | `src/solver/mod.rs:17-22` | — |
| Info | Public solver primitives rely on callers for nonzero-`w` guards | `src/solver` | — |

- **`AlmostSolved` budget (Low):** An `AlmostSolved` λ can be slightly dual-infeasible, yielding a slightly off budget `c = w·λ` with no post-solve re-check; tests note clarabel can land ~5e-4 off when the budget binds. Acceptable since the active-set loop re-validates `g` over the full grid each iteration; optionally verify `max_contact` before returning (helper exists at `refine_socp.rs:119-130`).
- **Debug-only consistency check (Low):** The only cross-check that recovered primal fuel matches the dual budget (5% band) is compiled out in release. Acceptable since residual is the load-bearing reachability check; if stronger guarantees are wanted, promote to a runtime `SolverFailed`/`NotConverged` on a looser tolerance.

### The three algorithms + orchestration

A faithful encoding of spec §5.1 with a clean `solve` / `solve_from_initial_times` split. No reachable panic or `unwrap` on legitimate input was found; the build and full suite (115 tests) pass. Findings are minor reporting/consistency nits.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | Reported residual measured on unpruned solution, not the pruned maneuver set returned | `src/algorithm/extract.rs:50-79` | — |
| Info | `solve_from_initial_times` rejects `n_init`/`n_coarse == 0` though that path never uses them | `src/algorithm/mod.rs:69-73,145-162` | — |
| Info | `refine`'s empty-candidate-set guard is effectively unreachable from the public API | `src/algorithm/refine.rs:80-84` | — |
| Info | Convergence target recomputed as `1+eps_cost` in two places (no drift; confirmed deliberate) | `src/algorithm/refine.rs:71-72,116-121` | — |

- **Residual vs pruned set (Low):** `Solution.residual` (line 56) is computed from the full min-fuel solution including maneuvers later pruned below `PRUNE_REL*max_dv`, so it describes a set that is not the one returned. The discrepancy is bounded (`<=1e-3` of the largest dv) and the split is deliberate/documented, but a consumer recomputing the residual of `Solution.maneuvers` can get a larger value. Either compute residual from the kept maneuvers, or add a one-line doc note on `types.rs Solution.residual`.

### Test suite & validation coverage

Strong, honest, and unusually well-reasoned; 115 tests green. The FD guards are the crown jewels, and the reported-not-asserted stance is scientifically defensible. Weaknesses are coverage gaps on error/failure paths — missing negative tests and one dead error variant, not correctness defects.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | `NotConverged` never exercised through the public `solve` API | `src/algorithm/mod.rs:121-162`; `tests/algorithm.rs` | — |
| Low | `KeplerDivergence` is a dead, unconstructible variant; Kepler loop can't report non-convergence | `src/types.rs:160-167`; `src/dynamics/kepler.rs:18-29` | — |
| Low | `SolverFailed` never tested at the public `solve()` level | `tests/algorithm.rs`; `tests/api.rs`; `src/solver/mod.rs:27-34` | — |
| Info | Monte-Carlo CI invariant 4 (Fig.8 shape) is a single weak inequality | `tests/monte_carlo.rs:118-124` | — |
| Info | `extract()` budget-vs-objective check is debug-only (untested in release) | `src/algorithm/extract.rs:44-48` | — |
| Info | `examples/mdot.rs` assertions never run under `cargo test` | `examples/mdot.rs:124-145` | — |

- **`NotConverged` untested at public API (Low):** Only triggered by the internal `refine()` unit test with `max_iters=1`; no integration test drives `solve()`/`solve_from_initial_times()` to a `NotConverged` result, so propagation through `run_pipeline` is untested. Add a test (or `Dynamics` mock) that cannot reach tolerance and assert `matches!(err, PlannerError::NotConverged { .. })`.
- **Dead `KeplerDivergence` (Low):** Declared but constructed nowhere; `mean_to_eccentric` returns the last iterate unconditionally. Either remove the variant or have the solver detect non-convergence and surface it; add a Kepler test at extreme `e` to document the convergence envelope. (Same root cause as the dynamics and cross-cutting findings — fix once.)
- **`SolverFailed` untested at public level (Low):** The wrapper-level mapping is covered, but no test drives the public pipeline into `SolverFailed` (e.g. a rank-deficient Γ on `T^opt`). Add one end-to-end test; low priority since the mapping itself is covered.

### Scientific / spec fidelity

A genuinely faithful, carefully validated encoding: every equation in spec §5 was cross-checked term-by-term and transcribed correctly, including the subtle points the spec flags. The honest disclosure of the paper's non-reproducible figures (author-confirmed typos) is sound. The one real gap is documentation drift in the normative math-reference section.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | Spec §5 math-reference describes the superseded Algorithm-3 QP, not the min-fuel SOCP the code runs | spec §5.1 (lines 210-220), §4.1 line 114, §4.2 line 162 | — |
| Info | "Faithful = match the paper's published numbers" claim in §1 contradicted by validated behavior | spec lines 31-33, 758-766; `Cargo.toml:8` | — |

- **Stale §5 Algorithm-3 (Low):** §5 is headed "Math reference (what the code must encode)" but its Algorithm 3 block specifies the paper's fixed-support-direction magnitude QP, while `solve()` → `extract()` (`src/algorithm/extract.rs:43`) calls `min_fuel_socp` — the Phase-5b direct min-fuel SOCP recovering full 3-DOF Δvⱼ, with the budget only a `debug_assert`. The deviation is documented, but only in the §6 changelog; a reader treating §5 as normative would encode the wrong algorithm. Update §5.1 (and the §4.1/§4.2 references) to describe the min-fuel SOCP `solve()` actually runs, presenting the printed QP as the paper's retained-primitive original, and note the change is a deliberate strong-duality-justified deviation.

### Monte Carlo validation harness

A high-quality, faithful reproduction of Fig.8 (three distinct seedings) and Fig.9 (solve time vs |T|). RNG is portable/deterministic, feature-gating is correct, and the CI invariants are principled and not over-fitted. No correctness defects; only minor robustness/style observations.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | Fig.8 summary/CDF grouping keys on the `n_init` count, would silently merge schemes sharing a count | `src/bin/monte_carlo.rs:256-277,344-354` | — |
| Info | Sampler logic duplicated between harness and CI test | `src/bin/monte_carlo.rs:122-134`; `tests/monte_carlo.rs:33-45` | — |
| Info | Fig.9 reports a single timed run per size (noisy point estimates) | `src/bin/monte_carlo.rs:375-409` | — |

- **Grouping by count (Low):** Correct today only because the three counts (2, 6, 10) are distinct; a future scheme reusing a count would be silently merged into one column. The scheme name is already stored on each row — group/summarize by the scheme label, or `debug_assert` the counts are unique.

### Project hygiene, packaging, CI, dependencies, docs structure

Well-organized for a research re-implementation: clean module tree, coherent `validation` feature flag with `dep:` syntax, dual MIT/Apache declared, CI covering fmt/clippy/build/test on a pinned toolchain, and healthy current dependencies with a committed lockfile. The gaps are documentation/packaging-oriented — the missing README and license texts are real friction for a crate intended to be cited and reused.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Medium | No README for a crate meant to be faithful, cited, and reused | repo root | — |
| Medium | SPDX declares `MIT OR Apache-2.0` but no license text files present | `Cargo.toml:6`; repo root | — |
| Low | CI does not enforce the declared MSRV (`rust-version = 1.92`) | `.github/workflows/ci.yml:17`; `Cargo.toml:5` | — |
| Low | rustdoc emits 29 unresolved intra-doc-link warnings that ship silently | `src/bin/monte_carlo.rs:68` and unit-annotation doc comments | — |
| Info | Transitive duplicate `thiserror` v1 (via clarabel) and a v1-era `bitflags` chain | `Cargo.lock` | — |
| Info | `version = 0.0.0` makes the crate non-publishable as-is | `Cargo.toml:3` | — |

- **Missing README (Medium):** No entry point for a reader/citer: no build instructions, no statement of the Φ₂₁ correction or the documented divergence from the paper's figures, no usage example. Add a README covering what/why, paper citation, the `validation` feature and how to reproduce Fig.8/9, the known divergences, and a minimal `solve` snippet.
- **Missing license texts (Medium):** The dual-license convention requires both full texts to be distributed; the SPDX field alone is not a license grant, so reusers/citers have no actual license to rely on. Add standard `LICENSE-MIT` and `LICENSE-APACHE` at the repo root — a 2-minute fix.
- **MSRV unenforced (Low):** `rust-version = "1.92"` is pinned but CI runs only on `@stable`, so an accidental newer-than-1.92 feature would pass CI yet break MSRV-honoring consumers. Add a `1.92` build job, or drop the pin if MSRV is not a commitment.
- **rustdoc warnings (Low):** Unit annotations like `[s]`/`[m]`/`[rad]` parse as intra-doc links, producing ~29 warnings that would render as broken-looking docs on docs.rs. Escape brackets or use backticks, and add `cargo doc --no-deps --all-features` to CI.

### Cross-cutting code quality

Unusually disciplined: zero panic/unsafe/todo in non-test code, no production `unwrap`/`expect` (two infallible documented ones aside), `Result`-based boundaries, centralized validation, and numerical hazards handled by formulation. The soft spots overlap with findings already raised in other dimensions.

| Severity | Title | Location | Status |
|----------|-------|----------|--------|
| Low | `KeplerDivergence` variant defined but unreachable — solver can't signal non-convergence | `src/types.rs:160-167`; `src/dynamics/kepler.rs:18-29` | — |
| Low | Unguarded GVE singularities in `B(t)` at degenerate inclinations | `src/dynamics/b_matrix.rs:20,29,32` | — |
| Low | Public `TimeGrid` API performs unvalidated `f64→usize` cast (saturation off the solve path) | `src/types.rs:50-67` | — |
| Info | `extract_qp` public but superseded by `min_fuel_socp`; only exercised via integration test | `src/solver/extract_qp.rs`; `src/lib.rs:16`; `tests/solver.rs:45` | — |

- **`TimeGrid` cast (Low):** `len()` computes `((t_f - t_i)/dt).round() as usize + 1` with no internal validation; the `as` cast saturates (negative ratio → 0; overflow → `usize::MAX`, after which `+1` wraps and `times()` attempts a catastrophic iteration). `TimeGrid::uniform`/`len`/`times` are public and callable directly without going through `validate_inputs`. Either document that these assume `dt>0 && t_f>=t_i`, or have `uniform` return a `Result`/`debug_assert` so misuse is caught rather than silently saturating.
- (`KeplerDivergence` and the `B(t)` `tan i` singularity are the same defects raised under Dynamics/Test suite — fix once and they clear in all three dimensions.)

## Prioritized recommendations

**P0 — must-fix correctness:** none. The audit found no Critical/High findings and no correctness defects in the shipped solve path.

> **Status (2026-06-19):** all five P1 items resolved — see
> `docs/superpowers/plans/2026-06-19-koenig-planner-p1-audit-fixes.md`. The P0
> set was empty; P2 items remain open.

**P1 — should-fix (fidelity, legal, and silent-misuse hardening):**

1. **Update spec §5 to describe the min-fuel SOCP** ("Spec §5 math-reference describes the superseded Algorithm-3 QP"). Make the normative math-reference match what `solve()` runs; present the printed QP as the retained paper-original, and qualify the §1 "faithful = match the numbers" claim per the Phase-5 reframing ("Faithful = match the published numbers" claim contradicted by validated behavior).
2. **Add `LICENSE-MIT` and `LICENSE-APACHE` files** ("SPDX declares MIT OR Apache-2.0 but no license text files present"). Resolves a real legal-grant gap.
3. **Add a README** ("No README for a crate meant to be faithful, cited, and reused"). State the Φ₂₁ correction, the documented paper-figure divergence, the `validation` feature, and a `solve` snippet.
4. **Resolve the `KeplerDivergence` dead variant and unguarded domains** (one fix clears three findings: "KeplerDivergence is a dead variant", "Kepler solver silently returns wrong value for e≥1", "Unguarded GVE singularities in B(t)"). Either make `mean_to_eccentric` return `Result` and emit `KeplerDivergence` on non-convergence, or remove the variant; in either case add `debug_assert`s for `0 ≤ e < 1` and inclination bounded away from `0`/`π`/`π/2`, or document those preconditions on `AbsoluteOrbit`/`J2Roe`.
5. **Guard or document the public `TimeGrid` cast and `Piecewise` phase assumption** ("Public TimeGrid API performs unvalidated f64→usize cast", "Piecewise silently assumes chief-at-apogee phase and t_i=0"). Add `debug_assert`/`Result` or prominent rustdoc preconditions so direct callers of these public types fail fast rather than silently producing garbage.

**P2 — nice-to-have / polish:**

6. **Clarify the residual contract** ("Reported residual measured on the unpruned solution") — compute from the kept maneuvers or add a doc note on `Solution.residual`.
7. **Tidy the superseded public surface** ("`extract_qp` still pub and re-exported", "`SolveParams.q` is a self-documented-dead field", "Crate-root re-exports expose implementation detail") — deprecate/relocate `extract_qp`, remove or mark `SolveParams.q`, and group core vs encoding types at the crate root.
8. **Strengthen negative-path tests** ("NotConverged never exercised through public API", "SolverFailed never tested at public level") — add end-to-end tests for both failure variants.
9. **Close CI/doc gaps** ("CI does not enforce the declared MSRV", "rustdoc emits ~29 intra-doc-link warnings") — add a `1.92` build job and a `cargo doc` step; escape unit-bracket annotations.
10. **Micro-optimizations & robustness** ("`vertex_columns()` recomputed per call" → hoist to `const`/`LazyLock`; "Fig.8 grouping keys on n_init count" → key on scheme label; "wrappers use clarabel default tolerances" → revisit only if degenerate grids start failing).
11. **Packaging when ready to share** ("`version = 0.0.0` makes the crate non-publishable") — bump to `0.1.0` and add `repository`/`readme`/`keywords`/`categories`; optionally add `cargo audit`/`cargo deny` to CI given the solver-heavy transitive tree.
12. **Update spec §4.1** to list `fuel_generator` on the `SublevelSet` sketch so the architecture reference matches the realized trait.

## Methodology & limitations

This audit covered ten dimensions: (1) architecture & module design, (2) J2 mean-ROE dynamics, (3) cost models, (4) convex-solver wrappers, (5) the three algorithms + orchestration, (6) test suite & validation coverage, (7) scientific/spec fidelity, (8) the Monte-Carlo harness, (9) project hygiene/packaging/CI/deps, and (10) cross-cutting code quality. Each dimension was a full read of the assigned source with cross-references and, where useful, building and running the crate. High-severity findings were slated for independent verification against the source; in this audit no Critical/High findings were raised, so the post-verification severity distribution is unchanged from the dimension reports (confirmed 0, partial 0, refuted 0). Equation-level fidelity claims were checked term-by-term against spec §5 and corroborated by the crate's own finite-difference oracles.

What this audit did **not** do: it did not re-derive Koenig & D'Amico's equations from the original IEEE TAC paper from scratch (fidelity was assessed against the design spec's §5 math reference, which prior phases verified against the paper and the author confirmed by DM); it did not run a code-coverage tool (coverage was assessed by reading the tests, not measured); it did not run `cargo audit`/`cargo deny` (not installed in the environment), so transitive-dependency advisory status is unverified; it did not execute the long Fig.9 sweeps to micro-benchmark performance; and it did not attempt to reproduce the paper's published worked-example figures (which prior phases established are not bit-reproducible due to author-confirmed typos in the source). The build/lint/test gate reported above was run live and observed (not fabricated); every file/line citation and numeric claim in this report was spot-checked against the source before publication.