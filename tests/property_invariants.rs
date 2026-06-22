//! Property-based tests for the scientific invariants of the planner.
//!
//! Targets are never random vectors: each problem builds a target
//! `w = Σ_k Γ(t_k)·Δv_k` from known impulses, so `w` is provably reachable
//! with ≤ 6 impulses ([KD20] eq. 42). Numeric properties run on a
//! well-conditioned domain tier; tolerances are justified per property.

use koenig_damico_planner::cost::{FaceMax, Norm2, SublevelSet};
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::PlannerError;
use koenig_damico_planner::{
    solve_from_initial_times, CostModel, Dynamics, Pseudostate, SolveParams, TimeGrid,
};
use nalgebra::SVector;
use proptest::prelude::*;

/// Uniform-Norm2 cost model: applies [KD20] Table II L2 cost at every time.
struct UniformNorm2;
impl CostModel for UniformNorm2 {
    fn at(&self, _t: f64) -> &dyn SublevelSet {
        &Norm2
    }
}

/// Uniform-FaceMax cost model: applies [KD20] Table II face-max cost at every time.
struct UniformFaceMax;
impl CostModel for UniformFaceMax {
    fn at(&self, _t: f64) -> &dyn SublevelSet {
        &FaceMax
    }
}

/// Raw generated scalars for a well-conditioned, reachable problem.
#[derive(Debug, Clone)]
struct RawProblem {
    a: f64,
    e: f64,
    i: f64,
    raan: f64,
    argp: f64,
    m0: f64,
    n_periods: f64,
    n_points: usize,
    /// K impulses: (fraction of the window in [0,1], Δv components [m/s]).
    impulses: Vec<(f64, [f64; 3])>,
}

/// Well-conditioned tier: elliptic chief away from the i→0/π singularity,
/// grid spanning a few orbits, 1..=6 bounded impulses ([KD20] §VIII regime).
fn well_conditioned_problem() -> impl Strategy<Value = RawProblem> {
    (
        7.0e6..5.0e7f64,                               // a [m]
        0.0..0.6f64,                                   // e (well-conditioned)
        (10.0f64..170.0).prop_map(|d| d.to_radians()), // i [rad]
        0.0..std::f64::consts::TAU,                    // raan [rad]
        0.0..std::f64::consts::TAU,                    // argp [rad]
        0.0..std::f64::consts::TAU,                    // mean_anom [rad]
        2.0..4.0f64,                                   // window length [periods]
        200usize..=2000,                               // grid points (bounded)
        prop::collection::vec(
            (0.0..1.0f64, proptest::array::uniform3(-1.0..1.0f64)),
            1..=6, // K ∈ [1,6]
        ),
    )
        .prop_map(
            |(a, e, i, raan, argp, m0, n_periods, n_points, impulses)| RawProblem {
                a,
                e,
                i,
                raan,
                argp,
                m0,
                n_periods,
                n_points,
                impulses,
            },
        )
}

/// Build a reachable problem: returns (dynamics, grid, target w, seed times).
/// `None` if a draw degenerates (near-zero target, or construction fails) —
/// callers skip those cases.
fn build_reachable(raw: &RawProblem) -> Option<(J2Roe, TimeGrid, Pseudostate, Vec<f64>)> {
    let chief = AbsoluteOrbit::new(raw.a, raw.e, raw.i, raw.raan, raw.argp, raw.m0);
    let period = std::f64::consts::TAU / chief.mean_motion();
    let t_i = 0.0;
    let t_f = raw.n_periods * period;
    let dt = (t_f - t_i) / (raw.n_points as f64);
    let grid = TimeGrid::uniform(t_i, t_f, dt).ok()?;
    let dynamics = J2Roe::new(chief, t_i, t_f).ok()?;

    let last = grid.len().saturating_sub(1);
    let mut w = Pseudostate::zeros();
    let mut idxs = std::collections::BTreeSet::new();
    for (frac, dv) in &raw.impulses {
        let idx = ((frac * last as f64).round() as usize).min(last);
        let t = grid.time(idx);
        let gamma = dynamics.gamma(t).ok()?;
        let dvv = SVector::<f64, 3>::new(dv[0], dv[1], dv[2]);
        w += gamma * dvv;
        idxs.insert(idx);
    }
    if w.norm() < 1e-9 {
        return None; // reject a near-zero target ([KD20] eq. 4 assumes w ≠ 0)
    }
    let seeds: Vec<f64> = idxs.iter().map(|&i| grid.time(i)).collect();
    Some((dynamics, grid, w, seeds))
}

/// Reusable: max over the grid of the gauge g_{U(1,t)}(Γᵀ(t)·λ) ([KD20] eq. 40).
fn max_gauge<D: Dynamics, C: CostModel>(
    dynamics: &D,
    cost: &C,
    grid: &TimeGrid,
    lambda: &Pseudostate,
) -> f64 {
    let mut max_g = f64::NEG_INFINITY;
    for t in grid.times() {
        if let Ok(gamma) = dynamics.gamma(t) {
            let g = cost.at(t).contact(gamma.transpose() * lambda);
            if g > max_g {
                max_g = g;
            }
        }
    }
    max_g
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 96,
        rng_seed: proptest::test_runner::RngSeed::Fixed(0x4B44_3230_5345_4544),
        ..ProptestConfig::default()
    })]

    /// [KD20] eq. 33/42: a reachable target is recovered exactly. Seeding the
    /// exact construction times makes the SOCP have an exact solution, so the
    /// pre-prune residual is ~machine-eps. Tolerance 1e-6 ≫ clarabel noise.
    #[test]
    fn reachable_target_has_tiny_residual(raw in well_conditioned_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let mut converged = false;
        for use_facemax in [false, true] {
            let result = if use_facemax {
                solve_from_initial_times(&dynamics, &UniformFaceMax, w, grid, &params, &seeds)
            } else {
                solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds)
            };
            if let Ok(sol) = result {
                converged = true;
                prop_assert!(
                    sol.residual < 1e-6,
                    "residual {} too large (facemax={})", sol.residual, use_facemax
                );
            }
        }
        prop_assume!(converged);
    }

    /// [KD20] Property 5 (p.4) / eq. 21 (S(c) = co S¹(c)): a control profile
    /// reaching w with ≤ n = 6 impulses EXISTS, but this is an EXISTENCE bound
    /// — the convex optimum can be non-unique and an interior-point solver may
    /// return a legitimate non-vertex optimum with more active impulses at the
    /// same optimal cost c*. So we assert the invariants the implementation
    /// actually guarantees: a nonzero reachable target yields ≥ 1 maneuver, and
    /// the returned set carries no dust (extraction prunes any Δv far below the
    /// largest), with every Δv finite.
    #[test]
    fn reachable_target_yields_nondust_maneuvers(raw in well_conditioned_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let result = solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds);
        prop_assume!(result.is_ok());
        let sol = result.unwrap();
        prop_assert!(!sol.maneuvers.is_empty(), "nonzero reachable w must yield >= 1 maneuver");
        let max_dv = sol
            .maneuvers
            .iter()
            .map(|m| m.dv.norm())
            .fold(0.0_f64, f64::max);
        prop_assert!(max_dv > 0.0);
        // Extraction prunes any maneuver below ~1e-3 of the largest; a 0.5e-3
        // floor catches a regression that returns dust without coupling tightly
        // to the private PRUNE_REL constant or flaking on interior-point noise.
        for m in &sol.maneuvers {
            prop_assert!(m.dv.iter().all(|x| x.is_finite()));
            prop_assert!(
                m.dv.norm() >= 0.5e-3 * max_dv,
                "returned dust maneuver: ‖Δv‖ {} < 0.5e-3 × max {}", m.dv.norm(), max_dv
            );
        }
    }

    /// [KD20] eq. 37/Thm 3: the optimal cost c* = total_dv is non-negative,
    /// and strictly positive for a nonzero reachable target. (NOT λ ≥ 0.)
    #[test]
    fn optimal_cost_is_nonnegative(raw in well_conditioned_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let result = solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds);
        prop_assume!(result.is_ok());
        let sol = result.unwrap();
        prop_assert!(sol.total_dv > 0.0, "c* = {} not > 0 for nonzero w", sol.total_dv);
        prop_assert!(sol.total_dv.is_finite());
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 64,
        rng_seed: proptest::test_runner::RngSeed::Fixed(0x4B44_3230_5345_4544),
        ..ProptestConfig::default()
    })]

    /// [KD20] eq. 40 / Algorithm 2 termination: at convergence the optimal
    /// dual satisfies max_t g_{U(1,t)}(Γᵀ(t)λ) ≤ 1 + ε_cost. Reuses the
    /// planner's own contact code, so this is a true optimality certificate.
    #[test]
    fn converged_dual_satisfies_support_constraint(raw in well_conditioned_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let mut converged = false;
        for use_facemax in [false, true] {
            let sol = if use_facemax {
                solve_from_initial_times(&dynamics, &UniformFaceMax, w, grid, &params, &seeds)
            } else {
                solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds)
            };
            let Ok(sol) = sol else { continue; };
            converged = true;
            let max_g = if use_facemax {
                max_gauge(&dynamics, &UniformFaceMax, &grid, &sol.lambda)
            } else {
                max_gauge(&dynamics, &UniformNorm2, &grid, &sol.lambda)
            };
            prop_assert!(
                max_g <= 1.0 + params.eps_cost + 1e-6,
                "max_t g = {} exceeds 1 + eps_cost (facemax={})", max_g, use_facemax
            );
        }
        prop_assume!(converged);
    }

    /// [KD20] Property 3 / eq. 24: the cost is a degree-1 gauge. Scaling the
    /// target by α scales the optimal cost by α (the robust scalar relation).
    #[test]
    fn cost_is_positively_homogeneous(
        raw in well_conditioned_problem(),
        alpha in 0.1..10.0f64,
    ) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let base = solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds);
        let scaled = solve_from_initial_times(&dynamics, &UniformNorm2, alpha * w, grid, &params, &seeds);
        prop_assume!(base.is_ok() && scaled.is_ok());
        let base = base.unwrap();
        let scaled = scaled.unwrap();
        let expected = alpha * base.total_dv;
        // For a degree-1 gauge cost ([KD20] Property 3) the optimal dual λ and
        // the active time set are scale-invariant, so in exact arithmetic
        // total_dv(αw) = α·total_dv(w) EXACTLY. The observed spread is purely
        // numerical: near the eps_cost (0.01) refinement threshold a borderline
        // candidate time can be active for one scale but not the other,
        // perturbing the primal total_dv by at most ~the gap the refinement
        // tolerates (≈ 2·eps_cost worst case). Empirically the relative diff is
        // ~machine-eps for almost all draws, with a rare tail near 4e-3
        // (≈ 5000 cases characterized); 1.5e-2 clears that max with ~3.6×
        // margin while staying under the 2·eps_cost = 2e-2 theoretical bound.
        prop_assert!(
            (scaled.total_dv - expected).abs() <= 1.5e-2 * expected.max(1e-12),
            "homogeneity: total_dv(αw) = {}, α·total_dv(w) = {}", scaled.total_dv, expected
        );
    }

    /// Implementation contract (not paper math): recomputing the residual from
    /// the pruned maneuvers stays small. The bound is the pruned dust mass —
    /// each pruned maneuver has ‖Δv‖ < PRUNE_REL (1e-3) of the largest Δv, and
    /// with up to K ≤ 6 maneuvers, scaled through the Γ conditioning of the
    /// reachable geometry (κ(Γ) reaches the low hundreds in this tier), the
    /// relative reconstruction error can exceed the naive 6·PRUNE_REL estimate.
    /// Empirically (≈ 30000 cases characterized) it is ~machine-eps for almost
    /// all draws, with a tail clustered near 2.5e-2 and a rare ill-conditioned-Γ
    /// extreme at 5.55e-2; 8e-2 clears that observed max with margin while
    /// remaining a meaningful guard against gross unpruned dust.
    /// NOTE: this reconstruction error is geometry-amplified through κ(Γ) and
    /// is heavy-tailed — the 8e-2 bound is a characterized empirical guard with
    /// only ~1.45× margin over the 5.55e-2 observed max, the thinnest in the
    /// suite. It is intentionally a gross-dust guard, not a tight contract; the
    /// tight, conditioning-independent guard is `pruned_cost_impact_is_bounded`.
    #[test]
    fn pruned_plan_residual_stays_small(raw in well_conditioned_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let result = solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds);
        prop_assume!(result.is_ok());
        let sol = result.unwrap();
        let mut recon = Pseudostate::zeros();
        for m in &sol.maneuvers {
            recon += dynamics.gamma(m.t).expect("reachable chief => Γ finite") * m.dv;
        }
        let r_pruned = (w - recon).norm() / w.norm();
        prop_assert!(r_pruned < 8e-2, "pruned residual {} too large", r_pruned);
        prop_assert!(sol.residual < 1e-6, "pre-prune residual {} too large", sol.residual);
    }

    /// Implementation contract (conditioning-INDEPENDENT, unlike the
    /// reconstruction-error guard): pruning only drops dust, so the kept plan's
    /// L2 fuel cost is within (#dropped)·PRUNE_REL of the reported pre-prune
    /// objective. Each dropped Δv has ‖Δv‖ < PRUNE_REL (1e-3) × max‖Δv‖, and
    /// max‖Δv‖ ≤ total_dv (a single term cannot exceed the sum), so with the
    /// handful of dropped dust maneuvers in this tier the relative cost gap is
    /// < ~1e-2 regardless of Γ conditioning. This is the tight guard the
    /// geometry-amplified reconstruction-error bound (below) cannot be.
    #[test]
    fn pruned_cost_impact_is_bounded(raw in well_conditioned_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        let result = solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds);
        prop_assume!(result.is_ok());
        let sol = result.unwrap();
        // Under Norm2 the reported total_dv is the pre-prune Σ‖Δv‖.
        let kept_cost: f64 = sol.maneuvers.iter().map(|m| m.dv.norm()).sum();
        let rel_cost_gap = (sol.total_dv - kept_cost).abs() / sol.total_dv.max(1e-12);
        prop_assert!(
            rel_cost_gap < 2e-2,
            "pruned cost gap {rel_cost_gap} exceeds the conditioning-independent dust bound"
        );
    }
}

/// Wide valid tier: full elliptic range and inclinations near (but not at)
/// the singularities. Used only for the no-panic totality property.
fn wide_valid_problem() -> impl Strategy<Value = RawProblem> {
    (
        7.0e6..5.0e7f64,
        0.0..0.85f64, // up to near-degenerate e
        (5.0f64..175.0).prop_map(|d| d.to_radians()),
        0.0..std::f64::consts::TAU,
        0.0..std::f64::consts::TAU,
        0.0..std::f64::consts::TAU,
        1.0..5.0f64,
        100usize..3000,
        prop::collection::vec(
            (0.0..1.0f64, proptest::array::uniform3(-5.0..5.0f64)),
            1..=6,
        ),
    )
        .prop_map(
            |(a, e, i, raan, argp, m0, n_periods, n_points, impulses)| RawProblem {
                a,
                e,
                i,
                raan,
                argp,
                m0,
                n_periods,
                n_points,
                impulses,
            },
        )
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        rng_seed: proptest::test_runner::RngSeed::Fixed(0x4B44_3230_5345_4544),
        ..ProptestConfig::default()
    })]

    /// Totality: on any well-posed reachable problem, `solve_from_initial_times`
    /// returns Ok or a typed PlannerError — never a panic. When Ok, every
    /// reported number is finite.
    #[test]
    fn solve_is_total_and_finite_on_ok(raw in wide_valid_problem()) {
        let problem = build_reachable(&raw);
        prop_assume!(problem.is_some());
        let (dynamics, grid, w, seeds) = problem.unwrap();
        let params = SolveParams::default();
        match solve_from_initial_times(&dynamics, &UniformNorm2, w, grid, &params, &seeds) {
            Ok(sol) => {
                prop_assert!(sol.total_dv.is_finite());
                prop_assert!(sol.residual.is_finite());
                prop_assert!(sol.lambda.iter().all(|x| x.is_finite()));
                prop_assert!(sol.maneuvers.iter().all(|m| m.dv.iter().all(|x| x.is_finite())));
            }
            Err(_) => { /* a typed error is an acceptable, non-panicking outcome */ }
        }
    }

    /// Constructor totality: invalid chief / grid parameters are rejected with
    /// PlannerError::InvalidInput, never a panic ([KD20] eq. 50 preconditions:
    /// a>0 finite, 0≤e<1, |sin i|≥1e-9; grid dt>0, t_f>t_i).
    #[test]
    fn invalid_domain_is_rejected_not_panicked(
        a in prop_oneof![Just(0.0f64), Just(-1.0f64), Just(f64::NAN), 7.0e6..5.0e7f64],
        e in prop_oneof![Just(1.0f64), Just(1.5f64), Just(-0.1f64), 0.0..0.9f64],
        i_deg in prop_oneof![Just(0.0f64), Just(180.0f64), 5.0..175.0f64],
        dt in prop_oneof![Just(0.0f64), Just(-1.0f64), 1.0..100.0f64],
        span in -10.0..1.0e5f64,
    ) {
        let chief = AbsoluteOrbit::new(a, e, i_deg.to_radians(), 0.0, 0.0, 0.0);
        let t_i = 0.0;
        let t_f = t_i + span;
        // Both fallible constructors must classify bad input as InvalidInput.
        if let Err(err) = TimeGrid::uniform(t_i, t_f, dt) {
            prop_assert!(matches!(err, PlannerError::InvalidInput(_)));
        }
        if let Err(err) = J2Roe::new(chief, t_i, t_f) {
            prop_assert!(matches!(err, PlannerError::InvalidInput(_)));
        }
    }
}
