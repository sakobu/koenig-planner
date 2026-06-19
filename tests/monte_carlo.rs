//! Phase 6 CI invariant test: Monte Carlo behavior of the public solver API on the
//! worked-example problem, driving the paper's THREE Fig. 8 initialization schemes
//! (n=2 endpoints, n=6 largest-g, n=10 evenly-spaced; Koenig & D'Amico 2020 p.11).
//! Asserts paper-INDEPENDENT invariants (NOT the paper's 4.90/3.99/3.31 means — see
//! spec §6 Phase 6 validation stance). Runs only under the `validation` feature (for
//! rand/rand_distr); CI runs `--all-features`.
#![cfg(feature = "validation")]

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_planner::{
    solve, solve_from_initial_times, PlannerError, Pseudostate, Solution, SolveParams, TimeGrid,
};
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

/// Solve one sample under the paper's seeding for column `k` (0 = n=2 endpoints,
/// 1 = n=6 largest-g, 2 = n=10 evenly-spaced) — mirroring the harness's `solve_scheme`.
fn solve_column(
    dynamics: &J2Roe,
    cost: &Piecewise,
    w: Pseudostate,
    grid: TimeGrid,
    k: usize,
) -> Result<Solution, PlannerError> {
    let p = SolveParams::default();
    match k {
        0 => solve_from_initial_times(dynamics, cost, w, grid, &p, &[grid.t_i, grid.t_f]),
        1 => solve(
            dynamics,
            cost,
            w,
            grid,
            &SolveParams {
                n_init: 6,
                ..SolveParams::default()
            },
        ),
        _ => {
            let times: Vec<f64> = (0..10)
                .map(|j| grid.t_i + (j as f64) * (grid.t_f - grid.t_i) / 9.0)
                .collect();
            solve_from_initial_times(dynamics, cost, w, grid, &p, &times)
        }
    }
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
        let (mut sum, mut count) = (0usize, 0usize);
        for &w in &ws {
            match solve_column(&dynamics, &cost, w, grid, k) {
                Ok(sol) => {
                    sum += sol.iterations;
                    count += 1;
                    max_iters = max_iters.max(sol.iterations);
                    max_res = max_res.max(sol.residual);
                }
                Err(_) => failures += 1,
            }
        }
        assert!(count > 0, "scheme n_init={n_init}: no successful solves");
        means[k] = sum as f64 / count as f64;
    }

    // Invariant 1: every solve succeeds (Phase 5b robustness).
    assert_eq!(failures, 0, "{failures} solve(s) failed; expected 0");
    // Invariant 2: converges within the paper's stated 8-iteration bound.
    assert!(
        max_iters <= 8,
        "max iterations {max_iters} exceeds the paper's 8-iter bound"
    );
    // Invariant 3: residual under 0.01% (the min-fuel SOCP reconstructs w).
    assert!(max_res < 1e-4, "max residual {max_res:.3e} exceeds 0.01%");
    // Invariant 4: Fig. 8 shape — the worst-case endpoints seed (n=2) needs more
    // refinement iterations than the well-spread evenly-spaced seed (n=10).
    assert!(
        means[0] > means[2],
        "mean iters n=2 endpoints ({:.3}) should exceed n=10 evenly-spaced ({:.3})",
        means[0],
        means[2]
    );

    eprintln!(
        "observed mean iters: n_init=2 -> {:.2}, 6 -> {:.2}, 10 -> {:.2}  (paper 4.90/3.99/3.31)",
        means[0], means[1], means[2]
    );
}
