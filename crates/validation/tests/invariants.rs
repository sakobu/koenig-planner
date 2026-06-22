//! Seeded CI invariant test for the Monte-Carlo harness. Asserts only paper-independent
//! invariants (the paper's means are a reference, not a target) and that the sampler shared
//! with the figure harness is deterministic and prefix-consistent. Ref: [KD20] Fig. 8 / Table V.

use koenig_damico_planner_validation::{
    run_fig8, sample_pseudostates, worked_example_cost, worked_example_dynamics, FIG8_SCHEMES, SEED,
};

const N_SAMPLES: usize = 64; // CI runtime (64 × 3 schemes = 192 solves)

#[test]
fn sampler_is_deterministic_and_prefix_consistent() {
    // Same (n, seed) → identical samples on every run/platform (portable StdRng).
    assert_eq!(
        sample_pseudostates(N_SAMPLES, SEED),
        sample_pseudostates(N_SAMPLES, SEED),
    );
    // The smaller draw is an exact prefix of the larger — the test and the figure harness
    // (N=200) share one sampler and one seed, so they cannot drift.
    let full = sample_pseudostates(200, SEED);
    let small = sample_pseudostates(N_SAMPLES, SEED);
    assert_eq!(&full[..N_SAMPLES], &small[..]);
}

#[test]
fn monte_carlo_invariants_hold() {
    let dynamics = worked_example_dynamics();
    let cost = worked_example_cost();
    let ws = sample_pseudostates(N_SAMPLES, SEED);
    let (rows, failures) = run_fig8(&dynamics, &cost, &ws, &FIG8_SCHEMES);

    assert_eq!(failures, 0, "{failures} solve(s) failed; expected 0");

    let max_iters = rows.iter().map(|r| r.iterations).max().unwrap_or(0);
    assert!(
        max_iters <= 8,
        "max iterations {max_iters} exceeds the paper's 8-iter bound"
    );

    let max_res = rows.iter().map(|r| r.residual).fold(0.0_f64, f64::max);
    assert!(max_res < 1e-4, "max residual {max_res:.3e} exceeds 0.01%");

    let mean = |n_init: usize| {
        let v: Vec<f64> = rows
            .iter()
            .filter(|r| r.n_init == n_init)
            .map(|r| r.iterations as f64)
            .collect();
        assert!(!v.is_empty(), "scheme n_init={n_init} produced no rows");
        v.iter().sum::<f64>() / v.len() as f64
    };
    let (m2, m6, m10) = (mean(2), mean(6), mean(10));
    // Worst-case 2-time endpoints needs more refinement than the 10-time evenly-spaced seed.
    assert!(
        m2 > m10,
        "mean iters n=2 ({m2:.3}) should exceed n=10 ({m10:.3})"
    );
    eprintln!("observed mean iters: 2->{m2:.2}, 6->{m6:.2}, 10->{m10:.2}  (paper 4.90/3.99/3.31)");
}
