# Koenig Planner — Rust Reimplementation Design

- **Date:** 2026-06-17
- **Status:** In implementation. **Phases 0–3 complete & verified** (CI green; Phase 1 dynamics confirmed across 5 independent routes — see `docs/superpowers/phase1-dynamics-verification-report.md`; Phase 2 cost models merged via PR #1, squash `51ac590`; Phase 3 solver wrappers merged via **PR #3**, squash `cdaff18`, +21 tests → 65, CI green); Phases 4–7 pending. **Resume at Phase 4 (three algorithms + orchestration).** See §6 for per-phase status.
- **Repo:** `github.com/sakobu/koenig-planner` (private). CI = GitHub Actions (`fmt` + `clippy -D warnings` + `build` + `test`, all `--all-features`); the Linux runner installs `libfontconfig1-dev` for the `plotters` validation feature.
- **Plans:** `docs/superpowers/plans/2026-06-17-koenig-planner-phase0-scaffolding.md`, `…-phase1-dynamics.md`, `…-phase2-cost-models.md`, `…-phase3-solver-wrappers.md`.
- **Source paper:** A. W. Koenig and S. D'Amico, "Fast Algorithm for Fuel-Optimal Impulsive Control of Linear Systems with Time-Varying Cost," *IEEE Transactions on Automatic Control*, 2020. DOI 10.1109/TAC.2020.3027804. (`docs/Planner.pdf`)

## 1. Goal & scope

Build a **faithful Rust reimplementation** of the paper's fuel-optimal impulsive
control algorithm. Faithful means: reproduce the algorithm exactly and match the
paper's published numbers.

**In scope (definition of done):**

1. The full three-step algorithm (Initialization, Iterative Refinement, Control-Input Extraction).
2. The J₂-perturbed mean relative-orbital-element (ROE) dynamics model (Appendix).
3. The two cost models used in the validation: `‖u‖₂` and `max(Vᶠᵃᶜᵉu)`, plus the
   time-varying piecewise selector (eq. 49).
4. **Worked-example validation**: Table III inputs → the Table IV 3-maneuver solution
   (≈ 82.4 mm/s, ~3 iterations) and the Fig. 7 contact-function curve.
5. **Monte Carlo validation** for the *proposed* algorithm: iteration-count
   distributions (Fig. 8) and compute-time-vs-discretization (Fig. 9).

**Out of scope (YAGNI):**

- Reference algorithms used only for comparison (Gilbert's algorithm, direct
  optimization — Table V). We reproduce *our* algorithm's timing/iteration data, not
  the head-to-head.
- Cost functions in Table I/II beyond the two the validation needs (`‖u‖₁`,
  `|u₁|+√(u₂²+u₃²)`). The `SublevelSet` trait makes these a later drop-in.
- `no_std` / flight-hardware targets, fixed-point, or deterministic-allocation work.
- A general LTV dynamics framework beyond the one J₂ ROE model (the `Dynamics` trait
  leaves the door open).

## 2. Decisions (locked)

| Decision | Choice |
|---|---|
| Primary goal | Faithful research reimplementation — correctness & matching the paper |
| Convex/QP solver | Native-Rust crate `clarabel` (covers SOCP, LP **and** QP) — no hand-rolled solver |
| Validation depth | Worked example **+** Monte Carlo (proposed algorithm only) |
| Structure | Single `koenig-planner` library crate, trait-based modules (Approach A) |

## 3. Background — the algorithm in brief

The planner drives a linear time-variant system from an initial state to a target
state at a fixed final time `t_f`, minimizing a (possibly time-varying) norm-like fuel
cost, using a small set of impulses (Δv's).

Key quantities:

- **Dynamics:** `ẋ = A(t)x + B(t)u`, with state transition matrix `Φ(t,t_f)`.
- **Pseudostate:** `w = x(t_f) − Φ(t_i,t_f)x(t_i)` — the part of the target the
  control must produce. `Γ(t) = Φ(t,t_f)·B(t)` maps an impulse at time `t` into
  pseudostate space.
- **Impulsive profile:** `u(t) = Σⱼ δ(t−tⱼ)vⱼ`, so `w = Σⱼ Γ(tⱼ)vⱼ`.

The problem (eq. 4) — minimize `∫ f(u,τ)dτ` s.t. `w = ∫ Γ(τ)u(τ)dτ` — is reformulated
via reachable-set theory into the **semi-infinite convex program (eq. 40):**

```
maximize_λ  λᵀw
subject to  max_{t∈T}  g_{U(1,t)}( Γᵀ(t)λ )  ≤  1
```

where `g_{U(1,t)}` is the **contact function** of the unit sublevel set of the cost at
time `t`. The optimal objective `c* = λ*ᵀw` is the minimum fuel cost (Theorem 3); `λ*`
is the outward normal to the reachable set at `w`.

## 4. Architecture (Approach A)

Single library crate `koenig-planner`; trait seams mirror the paper's math
abstractions.

```
koenig-planner/
├── Cargo.toml
├── src/
│   ├── lib.rs
│   ├── types.rs        # Pseudostate, Maneuver{t, dv}, TimeGrid, SolveParams, errors
│   ├── dynamics/
│   │   ├── mod.rs      # `Dynamics` trait
│   │   └── j2_roe.rs   # mean-element propagation, Φ(t,t_f), B(t), Γ(t)
│   ├── cost/
│   │   ├── mod.rs      # `SublevelSet` + `CostModel` traits
│   │   ├── norm2.rs    # ‖u‖₂  (unit ball)
│   │   ├── facemax.rs  # max(Vᶠᵃᶜᵉu)  (tetrahedral thrusters)
│   │   └── piecewise.rs# time-varying selector (eq. 49)
│   ├── solver/
│   │   ├── mod.rs
│   │   ├── refine_socp.rs  # builds + solves eq. 40 over a candidate-time set
│   │   └── extract_qp.rs   # Algorithm 3 QP
│   └── algorithm/
│       ├── mod.rs      # solve(...) orchestration + params
│       ├── init.rs     # Algorithm 1
│       ├── refine.rs   # Algorithm 2
│       └── extract.rs  # Algorithm 3
├── examples/mdot.rs            # worked example (Table III → Table IV, Fig. 7 data)
└── src/bin/monte_carlo.rs      # Fig. 8 / Fig. 9 harness
```

### 4.1 Trait interfaces

```rust
// N = 6 (ROE state), M = 3 (RTN Δv components).
pub const N: usize = 6;
pub const M: usize = 3;

// Dynamics: the algorithm only ever needs Γ(t) = Φ(t,t_f)·B(t).
pub trait Dynamics {
    fn gamma(&self, t: f64) -> SMatrix<f64, N, M>;
}

// The unit sublevel set U(1,t) of the cost at a given time.
// y = Γᵀ(t)·λ ∈ ℝᴹ.
pub trait SublevelSet {
    fn contact(&self, y: SVector<f64, M>) -> f64;             // g(y) = max_{z∈U} y·z
    fn support(&self, y: SVector<f64, M>) -> SVector<f64, M>; // s(y) = argmax z
    fn cone_constraints(&self, gamma_t: &SMatrix<f64, N, M>) -> ConicRows; // for SOCP/LP build
}

// Time-varying cost = piecewise selection of sublevel sets (eq. 49).
pub trait CostModel {
    fn at(&self, t: f64) -> &dyn SublevelSet;
}
```

`contact`/`support` are everything Algorithms 1–3 consume from the cost.
`cone_constraints` is what the solver layer needs to encode `g_{U(1,t)}(Γᵀ(t)λ) ≤ 1`
as SOC rows (`‖u‖₂`) or linear rows (`max(Vᶠᵃᶜᵉu)`).

### 4.2 End-to-end data flow

```
inputs: œ_c(t_i), w, CostModel, TimeGrid T, params(n_init, T^d, ε_cost, ε_remove, Q)
  → cache Γ(t) over T via Dynamics::gamma
  → [Init]    λ_est ∥ w; pick N times of largest g over T^d           ⇒ T^est
  → [Refine]  loop: solve eq.40 SOCP on T^est ⇒ λ_est; drop times <1−ε_remove;
              add local maxima of g >1; until max_t g ≤ 1+ε_cost      ⇒ (T^opt, λ_opt)
  → [Extract] directions sⱼ = s_{U(1,tⱼ)}; QP for magnitudes αⱼ        ⇒ Vec<Maneuver>
  → outputs: maneuvers [{t, Δv}], total Δv, iteration count, residual ‖w_err‖/‖w‖
```

## 5. Math reference (what the code must encode)

> **Transcription note:** the equations below are transcribed from the PDF for
> orientation. Phase 1 and Phase 2 MUST verify every term character-by-character
> against `docs/Planner.pdf` before trusting them, since the dynamics STM/B-matrix in
> particular are long and error-prone.
>
> **Validation status (2026-06-17).** Every equation, table, and validation target
> below was checked character-by-character against `docs/Planner.pdf`, and the
> high-risk `B(t)` and `Φ` were independently triangulated against their primary
> sources: Chernick & D'Amico 2018 (`docs/asm_2018_paper_chernickdamico.pdf`, the
> `B(t)` source), Koenig–Guffanti–D'Amico 2017 (`docs/Koenig_Guffanti_Damico.pdf`,
> ref [27], the STM source), and Hunter & D'Amico 2025
> (`docs/hunter_2025_ieee_aerospace_paper_final_v2.pdf`, an independent reproduction
> with its own validation case). No math errors were found. The inline notes below
> record the ambiguities that were resolved.

### 5.1 The three algorithms

**Algorithm 1 — Initialization.** Inputs `T^d` (coarse time samples), `λ_est`,
`Γ(t)`, `n_init` (user-specified initial candidate count). Compute
`g_{U(1,t)}(Γᵀ(t)λ_est)` for each `t ∈ T^d`; return the `n_init` times with the largest
values as `T^est`.

> **Notation:** `N`/`M` (= 6/3) are the *state/control dimensions* in the code; the
> paper's "N" for the number of candidate times is written `n_init` here, and the
> orbit-index "N ∈ ℤ" in eq. 49 is written `k` — all three are disambiguated to avoid
> the symbol collision.

**Algorithm 2 — Iterative Refinement.** Inputs `T^est`, `λ_est`, `T`, `w`, `Γ(t)`,
`ε_cost`, `ε_remove`.

```
do
  λ_est ← argmax λᵀw  s.t.  max_{t∈T^est} g_{U(1,t)}(Γᵀ(t)λ) ≤ 1      // eq. 40 on T^est
  for t ∈ T^est: if g_{U(1,t)}(Γᵀ(t)λ_est) < 1 − ε_remove: remove t   // drop slack times
  for local maxima of g_{U(1,t)}(Γᵀ(t)λ_est) over T:
        if value > 1: add t to T^est                                 // add violated times
while max_{t∈T} g_{U(1,t)}(Γᵀ(t)λ_est) > 1 + ε_cost
T^opt ← T^est;  λ_opt ← λ_est
```

`max_{t∈T} g` decreases monotonically toward 1, guaranteeing convergence.

**Algorithm 3 — Control-Input Extraction.** Inputs `T^opt`, `λ_opt`, `Γ(t)`, `w`,
`Q` (PD weight, identity in the example).

```
for tⱼ ∈ T^opt:
    sⱼ ← s_{U(1,tⱼ)}(Γᵀ(tⱼ)λ_opt)        // optimal direction
    yⱼ ← Γ(tⱼ)·sⱼ                          // its pseudostate contribution
α ← argmin  w_errᵀ Q w_err,  w_err = w − Σⱼ αⱼyⱼ
       s.t.  αⱼ ≥ 0,  Σⱼ αⱼ ≤ λ_optᵀw      // QP
u_opt = Σⱼ αⱼ·sⱼ  applied at tⱼ            // Maneuver{ t: tⱼ, dv: αⱼ·sⱼ }
```

### 5.2 Cost models (Table II)

| Cost `f(u)` | Contact `g(y)` | Support `s(y)` | Solver rows |
|---|---|---|---|
| `‖u‖₂` | `‖y‖₂` | `y/‖y‖₂` | one SOC: `‖Γᵀ(t)λ‖₂ ≤ 1` |
| `max(Vᶠᵃᶜᵉu)` | `maxₖ yᵀwₖ`, `wₖ ∈ cols(W)` | `argmaxₖ` column | linear: `wₖᵀΓᵀ(t)λ ≤ 1 ∀k` |

For `max(Vᶠᵃᶜᵉu)`: `W = [0_{M×1}, Vᵛᵉʳᵗᵉˣ]` (origin + thruster directions). From eq.
47/48 (tetrahedral, fixed-attitude occulter):

```
Vᵛᵉʳᵗᵉˣ = [[ √(2/3), −√(2/3),  0,      0     ],
           [ 0,       0,       √(2/3), −√(2/3)],
           [−√(1/3), −√(1/3),  √(1/3),  √(1/3)]]

Vᶠᵃᶜᵉ = (1/3)·[[−√(2/3), 0,      √(1/3)],
               [ √(2/3), 0,      √(1/3)],
               [ 0,     −√(2/3), −√(1/3)],
               [ 0,      √(2/3), −√(1/3)]]
```

**Sanity properties to test:** `λᵀs(λ) = g(λ)` (eq. 23); positive homogeneity
`g(αy) = α·g(y)` for `α ≥ 0`.

### 5.3 Piecewise time-varying cost (eq. 49)

```
f(u,t) = max(Vᶠᵃᶜᵉu)   for t ∈ T₁   (attitude-constrained windows)
       = ‖u‖₂          for t ∈ T₂
T₁ = { t : |t − (k+0.5)·T_orbit| < 1 hr,  k ∈ ℤ }   (2-hr windows around perigee)
T₂ = complement of T₁
```

### 5.4 Dynamics — J₂ mean-ROE model (Appendix)

> **Input convention — mean elements in, mean elements out.** `Φ` is evaluated at the
> chief's **mean** Keplerian elements and propagates **mean** ROE; `B` is the GVE matrix
> at the chief's mean elements, taken as the mean-ROE change under the near-identity
> Brouwer Jacobian approximation. Callers holding *osculating* elements must apply a
> Brouwer (first-order J₂ short-period) osc→mean conversion **before** calling the
> planner, and convert any osculating propagation back to mean before comparing. That
> transform lives in the I/O / validation layer, not the core planner; its closed form
> is in none of the source papers (use Brouwer / Vallado / Schaub if osculating inputs
> are ever needed).

Constants: `μ = 3.986e14 m³/s²`, `R_E = 6.378e6 m`, `J₂ = 1.082e-3`.

Absolute mean orbit `œ = [a, e, i, Ω, ω, M]`. Secular rates (eq. 50): `a, e, i`
constant;

```
Ω̇ = −(3 J₂ R_E² √μ)/(2 a^{7/2} η⁴) · cos i
ω̇ =  (3 J₂ R_E² √μ)/(4 a^{7/2} η⁴) · (5cos²i − 1)
Ṁ =  n + (3 J₂ R_E² √μ)/(4 a^{7/2} η³) · (3cos²i − 1),   n = √(μ/a³),  η = √(1−e²)
```

Propagate: `œ(t) = œ(t_i) + (t − t_i)·œ̇`.

**Kepler solve `M → E → ν`** (needed for `ν` in `B(t)`). No source paper writes this
out — Koenig, Chernick, Hunter and ref [27] all defer to "Kepler's equation" — so it
comes from standard astrodynamics (e.g. Vallado, *Fundamentals of Astrodynamics*):

```
M = E − e·sin E     Newton: E₀ = M + e·sin M;  E ← E − (E − e·sin E − M)/(1 − e·cos E)
                    reduce M to [−π, π);  iterate to ~1e-12 rad
ν = atan2( √(1−e²)·sin E ,  cos E − e )      quadrant-safe;  √(1−e²) = η
```

Well-conditioned at `e = 0.7` (`1 − e·cos E ≥ 1 − e = 0.3`), ~5–8 Newton iterations.
Because these relations are *not* in any PDF, they are verified by a round-trip
self-consistency test (`M→E→ν→E→M` identity) and known `M→ν` pairs, **not** by a
PDF cross-check (see Phase 1).

ROE state (eq. 51, chief = telescope, deputy = occulter, `Δ` = deputy − chief):

```
x = [ δa, δλ, δe_x, δe_y, δi_x, δi_y ]
  = [ Δa/a_c;
      ΔM + η_c(Δω + ΔΩ cos i_c);
      e_d cos ω_d − e_c cos ω_c;
      e_d sin ω_d − e_c sin ω_c;
      Δi;
      ΔΩ sin i_c ]
```

**Control input matrix `B(t)`** (`√(a_c/μ)` scaling; columns = R,T,N thrust on deputy;
`θ = ω + ν`, `ν` = true anomaly from `M` via Kepler):

> **Frame = RTN (≡ RIC), not NTW.** Per the paper: `R` is along the position vector
> (radial), `N` along the orbital angular-momentum vector (cross-track), and `T`
> completes the right-handed triad (`T = N × R`) — i.e. the *transverse* / along-track
> direction, **perpendicular to radial, not velocity-aligned**. At `e = 0.7` that gap is
> large near perigee, so this matters. It is the same triad as RIC (`T` = in-track,
> `N` = cross-track); it is **not** the NTW (velocity-tangent) frame. The `B`/`Γ`
> columns and the `u_R, u_T, u_N` outputs (Table IV) are all in this frame.

```
B(t) = √(a_c/μ) · [[B₁₁ B₁₂ 0  ],
                   [B₂₁ 0   0  ],
                   [B₃₁ B₃₂ B₃₃],
                   [B₄₁ B₄₂ B₄₃],
                   [0   0   B₅₃],
                   [0   0   B₆₃]]

B₁₁ = (2/η)e sin ν                         B₁₂ = (2/η)(1+e cos ν)
B₂₁ = −2η²/(1+e cos ν)                      B₃₁ = η sin θ
B₃₂ = η[(2+e cos ν)cos θ + e cos ω]/(1+e cos ν)
B₃₃ = η e sin ω sin θ / [tan i (1+e cos ν)] B₄₁ = −η cos θ
B₄₂ = η[(2+e cos ν)sin θ + e sin ω]/(1+e cos ν)
B₄₃ = −η e cos ω sin θ / [tan i (1+e cos ν)]
B₅₃ = η cos θ/(1+e cos ν)                   B₆₃ = η sin θ/(1+e cos ν)
```

**State transition matrix `Φ(t,t_f)`** — 6×6, quasi-nonsingular ROE with modified row 2
(δλ). Substitutions:

```
Δt = t_f − t,   κ = 3 J₂ R_E² √μ / (4 a^{7/2} η⁴),   G = η⁻²
P = 3cos²i − 1,  Q = 5cos²i − 1,  S = sin 2i,  T_sub = sin²i
e_{x1}=e cos(ω(t)),   e_{y1}=e sin(ω(t))
e_{x2}=e cos(ω(t_f)), e_{y2}=e sin(ω(t_f))
```

Nonzero terms (VERIFY against PDF in Phase 1):

```
Φ₁₁=1
Φ₂₁=(−1.5 n Δt − 7κηP)Δt   Φ₂₂=1   Φ₂₃=7κe_{x1}PΔt/η   Φ₂₄=7κe_{y1}PΔt/η   Φ₂₅=−7κηSΔt
Φ₃₁=3.5κe_{y2}QΔt          Φ₃₃=cos(ω̇Δt)−4κe_{x1}e_{y2}GQΔt
Φ₃₄=−sin(ω̇Δt)−4κe_{y1}e_{y2}GQΔt                       Φ₃₅=5κe_{y2}SΔt
Φ₄₁=−3.5κe_{x2}QΔt         Φ₄₃=sin(ω̇Δt)+4κe_{x1}e_{x2}GQΔt
Φ₄₄=cos(ω̇Δt)+4κe_{y1}e_{x2}GQΔt                         Φ₄₅=−5κe_{x2}SΔt
Φ₅₅=1
Φ₆₁=3.5κSΔt   Φ₆₃=−4κe_{x1}GSΔt   Φ₆₄=−4κe_{y1}GSΔt   Φ₆₅=2κT_sub Δt   Φ₆₆=1
```

> **Φ₂₄ is intentionally nonzero** (`Φ₂₄ = 7κe_{y1}PΔt/η`): under J₂ the δλ row couples
> to δe_y. The printed 6×6 matrix box in `Planner.pdf` shows `0` at cell (2,4) — that is a
> rendering/transcription artifact; the paper's own term list (which this block follows)
> is authoritative. Confirmed nonzero by the STM primary source ref [27]
> (Koenig–Guffanti–D'Amico 2017), Chernick & D'Amico 2018 eq. (32), and Hunter & D'Amico
> 2025 eq. (76).
>
> **The δλ row-2 modification lives in `Φ` only, never `B`.** Koenig's
> δλ = ΔM + η_c(Δω + ΔΩ cos i_c) (eq. 51) differs from the standard quasi-nonsingular δλ
> of ref [27]; the modification is absorbed **entirely into Φ row 2** — the `/η` on
> Φ₂₃/Φ₂₄ and the modified Φ₂₁. `B(t)` needs **no** change: its δλ row stays
> `[−2η²/(1+e cos ν), 0, 0]`, because the only place the modification could surface is
> B₂₃ (the ΔΩ / out-of-plane coupling), which is 0 under the independent-out-of-plane-
> maneuver approximation. So `Γ = Φ·B` is self-consistent exactly as transcribed.
>
> **Coefficient provenance (verified directly against ref [27] eq. (A6), 2026-06-17).**
> Reproduce **Koenig 2020's** row 2 verbatim (`Φ₂₁=(−1.5nΔt−7κηP)Δt`, `Φ₂₃=7κe_{x1}PΔt/η`,
> `Φ₂₄=7κe_{y1}PΔt/η`, `Φ₂₅=−7κηSΔt`). This is what produced Table IV, so **any other
> row-2 form will miss the published 82.4 mm/s solution.** The STM primary source ref [27]
> eq. (A6) uses the *standard* δλ and so prints row 2 differently —
> `Φ₂₁=−(1.5n+3.5κ·E·P)Δt`, `Φ₂₃=κe_{xi}·F·G·PΔt`, `Φ₂₄=κe_{yi}·F·G·PΔt`, `Φ₂₅=−κ·F·SΔt`,
> with `E=1+η`, `F=4+3η`, `G=η⁻²` (ref [27] eqs. 13–14). Koenig collapses these into the
> `/η` form above when he swaps in his η-weighted δλ; the difference is the **documented**
> row-2 modification, not a transcription error. Rows 1, 3, 4, 5, 6 are identical between
> Koenig and ref [27]. (`F=4+3η` is the primary-source value; Hunter eq. (77)'s `4+2η` is a
> Hunter typo, irrelevant here since Koenig's modified row uses no `F`.)

Then `Γ(t) = Φ(t,t_f)·B(t)`.

### 5.5 Units & scaling convention

- **State & matrices stay native / dimensionless.** `x` (eq. 51) is dimensionless
  (`δa = Δa/a_c`); `Φ` is dimensionless; `B(t) = √(a_c/μ)·[Bᵢⱼ]` maps an impulse `Δv`
  [m/s] to a *dimensionless* ROE change. **`√(a_c/μ)` equals Chernick/Hunter's
  `1/(n_c·a_c)`** (since `n = √(μ/a³)`) — it is *not* an extra factor; do not
  "simplify" it away. Hence `Γ = Φ·B` also maps `Δv` [m/s] → a dimensionless pseudostate.
- **Apply `a_c` exactly once, at the I/O boundary.** The native pseudostate
  `w = x(t_f) − Φ(t_i,t_f)x(t_i)` is dimensionless. Table III lists `w` in metres only
  because it is pre-multiplied by `a_c` for display
  (`[50, 5000, 100, 100, 0, 400] m = a_c·w_nd`). On input, divide by `a_c` (= 25000 km)
  to get `w_nd` before solving; multiply states by `a_c` for reporting. Never bake `a_c`
  into `B` or `Φ`.
- **Dual `λ`.** `λ` lives in the dual of whichever `w`-units you choose. Only its
  *direction* and the scale-invariant ratio `λᵀw / g(λ)` matter, so the `w`-scaling
  choice does not change optimal times/directions/Δv. Koenig's `λ_opt ≈ 1e-6·[…]` pairs
  with the metre-scaled `w` (its ~1e-6 magnitude × metre-scale `w` → mm/s Δv).

## 6. Roadmap (phases & exit criteria)

### Phase 0 — Scaffolding ✅ Done (2026-06-17, commits `2c6dec1`…`1e43195`)
Cargo lib + deps; `types.rs` (`Pseudostate = SVector<f64,6>`, `Maneuver{t, dv:SVector<f64,3>}`,
`TimeGrid`, `SolveParams`, error enum); empty trait defs compile.
**Exit:** `cargo test` green on stubs; CI runs.
**Done:** crate rooted at the repo root; full module tree + trait seams (`Dynamics`,
`SublevelSet`, `CostModel`) + stubs compile; 12 tests green; GitHub Actions CI green.
Deps locked: `nalgebra 0.35`, `clarabel 0.11`, `thiserror 2.0`, `approx 0.5` (dev),
`csv 1.4`/`plotters 0.3` (optional, behind a `validation` feature).

### Phase 1 — J₂ mean-ROE dynamics (highest correctness risk) ✅ Done & verified (commits `f2c7277`…`55a14dd`)
Mean-element secular propagation (eq. 50); Kepler solve `M→E→ν` (Newton, must handle
`e = 0.7`); `B(t)`; ROE STM `Φ(t,t_f)`; `gamma(t) = Φ·B`.
**Tests:** `Φ → I` as `Δt → 0`; finite-difference cross-check of `B(t)` columns vs
numeric `∂x/∂(Δv)`; a few hand-computed reference values; dimensional sanity;
**Kepler round-trip `M→E→ν→E→M` identity + known `M→ν` pairs** (the Kepler relations are
not in any source PDF — verify by self-consistency, not PDF cross-check).
**Exit:** all dynamics tests pass; every STM/`B` term verified character-by-character
against `docs/Planner.pdf` — *except* the Kepler block, which has no paper counterpart
and is covered by the round-trip / known-value tests instead. The `Φ₂₄`-nonzero and
`a_c`-scaling conventions (§5.4–5.5) are locked before anything depends on Phase 1.
**Done:** files `dynamics/{constants,kepler,orbit,b_matrix,stm,j2_roe}.rs`; **27 tests green**
(incl. the independent FD test `tests/fd_b_matrix.rs`). **Verified across 5 independent routes,
zero discrepancies** (full evidence: `docs/superpowers/phase1-dynamics-verification-report.md`):
(1) entrywise vs a Python oracle transcribed from §5.4 (Rust ↔ Python ≤ 1e-9);
(2) **character-by-character read of every `B`/`Φ`/secular term vs `docs/Planner.pdf`** — the
named exit gate — done by hand *and* by a 10-agent adversarial pass;
(3) primary-source triangulation: `B` vs Chernick & D'Amico 2018 **Eq. (38)**, `Φ` vs
Koenig–Guffanti–D'Amico 2017 **A6** (rows 1/3/4/5/6 identical; row 2 = the documented `/η` δλ
modification, confirmed against ref [27]'s `E=1+η` form);
(4) independent **finite-difference** of `B(t)` via a Cartesian r,v route using none of `B`'s
formulas (≤ 1e-9);
(5) the locked-convention tests (`Φ→I`, `Φ₂₄`≠0, `√(a/μ)` law, Kepler round-trip, anchors).
Subtleties independently confirmed: `Ṁ` uses `η³` (not `η⁴`); `Φ`'s mixed initial/final
eccentricity subscripts (`Φ₃₃: e_{x1}e_{y2}`, `Φ₄₄: e_{y1}e_{x2}`); `B₃₃` positive vs `B₄₃`
negative. Bug caught + fixed: at `M=π`, `wrap_to_pi(π)=−π` ⇒ `ν=−π` (apoapsis) — Kepler test
compares `|ν|`.

### Phase 2 — Cost models (Table II, eq. 47–49) ✅ Done & verified (PR #1, squash `51ac590`)
`Norm2`, `FaceMax` (with `Vᵛᵉʳᵗᵉˣ/Vᶠᵃᶜᵉ`), `Piecewise` selector.
**Tests:** `λᵀs(λ)=g(λ)`; positive homogeneity; known directions; `T₁/T₂` window logic.
**Exit:** cost tests pass.
**Done:** files `cost/{norm2,facemax,piecewise}.rs` + `cost/mod.rs` wiring test; **18 cost tests green**
(9 FaceMax + 6 Norm2 + 2 Piecewise + 1 wiring; full suite 44, CI green). The `SublevelSet`/`CostModel`
traits are **fully implemented — no `unimplemented!` left in the cost layer**. Key choices:
(1) **`cone_constraints` implemented now, not deferred** — `Norm2` → one SOC row `(Γᵀ, 1)`; `FaceMax`
→ four linear rows `(Γvₖ, 1)`; the origin column of `W=[0|V_vertex]` is vacuous as a cone row and is
omitted. Phase 3's `refine_socp` consumes these rows directly. (2) **`Piecewise::new(period)` takes the
orbit period as an argument** — the paper leaves Keplerian-vs-perturbed open; the worked example passes
`2π/n` (≈10.93 hr, paper rounds to 10.92). (3) **`Piecewise` dropped `Default`** (now carries fields) →
`cost/mod.rs` wiring test constructs via `Piecewise::new`. (4) **Citation fix:** positive homogeneity is
**eq. 8 / Property 3**, not eq. 23 (eq. 23 = the support identity only). `V_face` appears only as a
test-module transcription cross-check (`f(vₖ)=1/9` for every vertex); the algorithm uses `V_vertex`/`W`
exclusively. Verified three independent ways: gate-tested in a throwaway worktree before the plan was
finalized; an adversarial agent re-applied the plan's code and re-ran the full CI gate (0 blockers) and
reproduced every reference number; and the in-repo gate ran per task + CI on the PR. Reference numbers
from an independent pure-Python oracle; the cost matrices (eq. 47–49, Table II, eq. 23) re-verified
character-by-character against `docs/Planner.pdf`. Two gate gotchas caught + fixed: `clippy::approx_constant`
on the `1/√2` literal (L2 test uses the `(3,4,12)→13` vector); the exact eq. 49 window boundary
(`|t−center|=3600 s`) is a floating-point knife-edge (the boundary test probes ±1 s either side).

### Phase 3 — Solver wrappers ✅ Done & verified (merged via PR #3, squash `cdaff18`)
`refine_socp`: assemble eq. 40 over a candidate-time set into clarabel conic form
(linear + SOC cones from each time's `cone_constraints`), map maximize→minimize, return
`λ` + objective. `extract_qp`: the Algorithm 3 QP.
**Tests:** small hand-checkable problems with closed-form optima.
**Exit:** solver tests pass.
**Done:** files `src/solver/{mod,refine_socp,extract_qp}.rs` + re-exports in `lib.rs`; integration
test `tests/solver.rs`. **+21 tests → 65 total, CI gate green** (`fmt` + `clippy --all-features -D warnings`
+ `build` + `test`; the original +17 plus 4 added by a post-merge self-audit: `w=0`, unbounded-SOCP →
`SolverFailed` through the wrapper, and negative/NaN-budget rejection). Both wrappers are stateless/pure — they consume pre-assembled data, so the
`Dynamics`/`CostModel` traits are not referenced here. **clarabel 0.11.1 API pinned** (from the installed
crate source): solves `min ½xᵀPx+qᵀx s.t. Ax+s=b, s∈K`; `DefaultSolver::new(&P,&q,&A,&b,&cones,settings)
→ Result`; read `solution.{x,status,obj_val}`; **`solver.solve()` requires `use clarabel::solver::IPSolver`
in scope** (clarabel's own examples pull it via a `use …::*` glob — an explicit import list MUST add it);
cones via `use clarabel::solver::*` (`NonnegativeConeT`/`SecondOrderConeT`); `CscMatrix::from(&dense)` drops
exact zeros; **`P` must be passed UPPER-TRIANGULAR** (clarabel never symmetrizes — `kkt_assembly.rs`:
"user provided P is always triu regardless"). **Encodings (independently re-derived + 4-agent-verified
before coding):** *SOCP* — `x=λ`, `P=0`, `q=−w` (maximize→minimize), recover `c*=w·λ` directly from the
primal (not `obj_val`); `Norm2` time → one `SecondOrderConeT(M+1=4)` block `A=[0ᵀ;−Γᵀ]`, `b=[1,0,0,0]`;
`FaceMax` time → four `NonnegativeConeT` rows `(Γvₖ)ᵀλ≤1`; rows assembled **linear-first then SOC** in
lockstep with the cone vector; **`λ` is FREE (no sign cone — do not copy the QP's `α≥0`)**. *QP* —
`P=2·YᵀQY` (triu), `q=−2·YᵀQw`, drop const `wᵀQw` (add back for residual), one `NonnegativeConeT(K+1)`
with `A=[−I_K;1ᵀ]`, `b=[0_K;budget]`, `budget=λ_optᵀw`; `Q` symmetrized defensively. **Status mapping:**
accept `Solved`+`AlmostSolved`; every other `SolverStatus` → `PlannerError::SolverFailed` naming the status.
**KEY DECISION:** `refine_socp` returns `{lambda, objective}` only — the **"per-time slack" is the caller's
job** (Phase 4 recomputes `g_{U(1,t)}(Γᵀ(t)λ)` via `SublevelSet::contact`; it scans `g` over the full grid
`T` anyway, and `refine_socp` consumes `ConicRows`, not cost objects). Closed-form optima reproduced to
`1e-6`: pure-SOC `c*=13`; **face-max LP `c*=√3`** (Risk 6 — was un-cross-validated by the literature, here
hand-derived); mixed SOC+LP `c*=√3+1` (validates cone ordering = the realistic eq.49 Piecewise case); QP
α-cases covering interior / budget-binding / nonneg-binding / weighted-`Q` / **non-orthogonal-`Y`** (guards
the triu-`P` packing + factor-of-2 — a diagonal-only test can't) / singular-`P` (assert the unique residual,
not the non-unique `α`). **Two plan gaps caught during execution (both fixed):** (1) the `IPSolver` import
above; (2) helpers used only by `#[cfg(test)]` trip `dead_code` under `clippy -D warnings` until their
consumer lands — use `#[allow(dead_code)]` transiently, **not** `#[expect]` (which mis-fires in the
`cfg(test)` build where the tests *do* use them). **Phase 4 hand-off:** build `Vec<ConicRows>` via
`cost.at(t).cone_constraints(&dynamics.gamma(t))` for each `t ∈ T^est` → `refine_socp(w, &rows)`; then
`sⱼ = cost.at(tⱼ).support(Γᵀ(tⱼ)λ)`, `yⱼ = gamma(tⱼ)·sⱼ`, `extract_qp(w, &ys, &Q, budget=λ·w)` → `αⱼ`,
emit `Maneuver{ t: tⱼ, dv: αⱼ·sⱼ }`; **filter zero-support times before extract** (a `yⱼ=0` column leaves
`αⱼ` irrelevant-but-unconstrained).

### Phase 4 — Three algorithms + orchestration
Alg. 1 init; Alg. 2 refine (incl. discrete local-maxima finder over the grid, with
grid-endpoint handling); Alg. 3 extract; `solve(...)` wiring with Γ(t) caching.
**Tests:** `max_t g` decreases monotonically across iterations; convergence within
`1+ε_cost`; small residual on a synthetic case.
**Exit:** end-to-end `solve` runs on a synthetic problem and converges.

### Phase 5 — Worked-example validation (`examples/mdot.rs`)
Encode Table III + params; reproduce the published result. **Targets** in §7.
**Exit:** all §7 worked-example assertions pass within stated bands.

### Phase 6 — Monte Carlo harness (`src/bin/monte_carlo.rs`)
200 pseudostates `~ N(0, σ=1km)` per ROE; record iterations + wall-time. Reproduce
Fig. 8 (iteration CDFs for 2/6/10-time inits) and Fig. 9 (compute time vs `|T|`,
10→10⁶). Emit CSV (+ optional `plotters` PNGs).
**Exit:** iteration means near 4.90 / 3.99 / 3.31; Fig. 9 shape (≈constant ≤10⁴, then
linear) reproduced; residuals `< 0.01%` across the 200 cases.

### Phase 7 — Polish
Rustdoc cross-referencing code ↔ equation numbers; README; runnable examples; CI green.

## 7. Validation targets (precise numbers)

**Worked example inputs (Table III):**
- Initial mean absolute orbit: `a=25000 km, e=0.7, i=40°, Ω=358°, ω=0°, M=180°`.
- Target pseudostate (scaled by `a`): `[aδa, aδλ, aδe_x, aδe_y, aδi_x, aδi_y] = [50, 5000, 100, 100, 0, 400]` m.
- `t_i = 0`, `t_f = 117990 s` (3 orbits; `T_orbit = 10.92 hr`).
- Control time domain `T`: uniform 30 s grid over `[t_i, t_f]` → **3934 candidate times**.
- Cost: eq. 49 (FaceMax in 2-hr perigee windows, Norm2 elsewhere).
- Params: `T^d = 20` evenly spaced times, `n_init = 6`, `λ_est ∥ w`, `ε_cost = ε_remove = 0.01`, `Q = I`.

**Worked example expected outputs:**
- **3 maneuvers** at `t = [16050, 23280, 107100] s`:
  - `u_R = [9.68, 0.00, 16.51]` mm/s
  - `u_T = [−23.02, −0.40, 15.68]` mm/s
  - `u_N = [−25.56, −0.04, 40.26]` mm/s
- Total Δv ≈ **82.4 mm/s** (≤ 1% above the 82.0 mm/s lower bound).
- `λ_opt ≈ 1e-6 · [34.97, 3.42, 30.68, 17.84, −9.34, 146.79]ᵀ`.
- **~3 iterations** of Algorithm 2.
- Residual `‖w_err‖₂/‖w‖₂ < 0.01%`.

**Monte Carlo targets:**
- 200 pseudostates from zero-mean Gaussian, σ = 1 km per ROE.
- Mean iteration counts: **4.90 / 3.99 / 3.31** for 2 / 6 / 10-time initialization.
- Convergence within 1% of optimum in ≤ 8 iterations across all cases.
- Compute time roughly constant for `|T| ≤ 10⁴`, then linear (Fig. 9).

**Second worked example — independent cross-check** (Hunter & D'Amico 2025, "Sequential
Formulation Validation", *identical* J₂ ROE dynamics; use as additional integration
coverage, the Table III/IV case stays primary):
- Chief mean orbit: `a = 25000 km, e_x = −0.658, e_y = −0.239` (`e ≈ 0.70, ω ≈ 200°`),
  `i = 51°, Ω = 30°`, `u₀ = 65°`.
- Control window `39000 s` (~1 orbit); uniform 10 s grid → **3901 candidate times**.
- Target pseudostate `w = [0.66, −1.52, −0.38, −1.44, 0.29, −0.91] m`.
- Koenig-solver expected output: **3 maneuvers, total Δv ≈ 23.03e-5 m/s, ~4 iterations**
  (`ε_cost = ε_remove = 0.01`); dual lower bound 22.94e-5 m/s; residual `< 0.01%`.

> Assertions use sensible numerical bands (Δv, residual, iteration count), **not
> bit-for-bit equality** — exact figures depend on solver tolerances. The paper's
> tolerances are the spec.

## 8. Dependencies

| Crate | Role |
|---|---|
| `nalgebra` | static-dim linear algebra (`SMatrix<6,3>`, `SVector`) |
| `clarabel` | native-Rust conic solver — SOCP/LP (refinement) **and** QP (extraction) |
| `thiserror` | error types |
| `approx` (dev) | float-tolerant test assertions |
| `csv`, `plotters` (validation bin only) | emit Fig. 7/8/9 data |

## 9. Risks & mitigations

1. **J₂ ROE STM / B-matrix transcription (Phase 1)** — long, error-prone. *Mitigate:*
   character-by-character verification against the PDF + the Phase 1 property tests;
   build and lock Phase 1 before anything depends on it.
2. **Solver convention mismatch** (maximize vs minimize, cone ordering in clarabel).
   *Mitigate:* Phase 3 closed-form unit tests before integration. **[✅ retired in Phase 3:**
   `q=−w` + SOC layout `[0ᵀ;−Γᵀ]`/`b=[1,0,0,0]` verified against clarabel's shipped `example_socp.rs`
   and reproduce all closed-form optima to `1e-6`; `λ` is free; rows linear-first then SOC.**]**
3. **Exact-number matching** depends on solver tolerances. *Mitigate:* assert on bands
   tied to the paper's tolerances, not bit-equality.
4. **Discrete local-maxima finder (Alg. 2)** edge cases at grid ends / plateaus.
   *Mitigate:* dedicated unit tests.
5. **Kepler convergence at `e = 0.7`.** *Mitigate:* Newton with a robust initial guess;
   test against known `M→ν` pairs. (The explicit `M→E→ν` relations are not in any source
   PDF — see §5.4; verify by round-trip self-consistency, not paper cross-check.)
6. **Face-max (LP) cost path has no published reference optimum.** Koenig and Hunter
   validate only the L2/SDP cost cases; the `max(Vᶠᵃᶜᵉu)` LP path is un-cross-validated by
   the literature. *Mitigate:* a standalone Phase-3 unit test with a hand-derived
   closed-form LP optimum (not only the eq. 49 piecewise case, where the LP and L2 rows
   are entangled). **[✅ retired in Phase 3:** `s2_face_max_lp_closed_form` asserts the standalone
   hand-derived optimum `c*=√3` (`w=(0,0,1,0,0,0)`, `Γ=[I₃;0]`), disentangled from L2.**]**
7. **Non-smooth contact × local-maxima finder.** For the face-max cost the contact
   `maxₖ yᵀwₖ` is evaluated by enumerating the columns of `W` (no smooth closed form),
   which interacts with risk #4. *Mitigate:* exercise the Alg. 2 grid local-maxima finder
   specifically on the face-max cost.
8. **`clarabel` is neither source paper's solver** (Koenig used MATLAB + CVX; Hunter used
   CVX + SDPT3). Both confirm the problem class — Hunter independently re-derives the same
   eq. 40 SOCP (its eqs. 18–20) and reuses Koenig's solver as a black box — but exact
   digits depend on `clarabel`'s interior-point tolerances. *Mitigate:* band assertions
   tied to `ε = 0.01` (risk #3); optionally warm-start `clarabel` across refinement
   iterations (the per-iteration SOCPs share `λ_opt`).

## 10. Open questions

None blocking. The dynamics ambiguities raised during validation are resolved and
documented inline (§5.4 `Φ₂₄` nonzero, δλ/`B` consistency, mean-element convention;
§5.5 scaling; §5.4 Kepler). Other Table I/II cost functions and a second dynamics model
are deferred behind the existing traits and can be added without architectural change.

*Remaining nice-to-haves (non-blocking):* a standard astrodynamics text (Vallado) is the
citable home for the Kepler relations; and the `4+Nη`-style packaging of the modified
`Φ` row-2 differs between Hunter and Chernick — irrelevant to a faithful Koenig
reproduction, which uses Koenig's verbatim `Φ` terms (confirmed against ref [27]).
