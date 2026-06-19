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
        let params = SolveParams {
            n_init,
            ..SolveParams::default()
        };
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
    assert!(
        max_iters <= 8,
        "max iterations {max_iters} exceeds the paper's 8-iter bound"
    );
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
