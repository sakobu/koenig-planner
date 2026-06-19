//! Monte Carlo harness — Fig. 8 (Algorithm-2 iteration distributions) and Fig. 9
//! (compute time vs discretization |T|) for the *proposed* algorithm on the
//! worked-example problem (Table III chief, eq. 49 cost).
//!
//! Build & run with the `validation` feature (needs rand/csv/plotters):
//!
//!   cargo run --features validation --bin monte_carlo            # both sweeps
//!   cargo run --features validation --bin monte_carlo -- fig8    # Fig. 8 only
//!   cargo run --features validation --bin monte_carlo -- fig9    # Fig. 9 only
//!
//! Validation stance (spec §6 Phase 6): we REPORT our iteration distributions and
//! compare to the paper's 4.90/3.99/3.31 means as *reference*, not bit-reproduce the
//! paper's (internally inconsistent) figures. See tests/monte_carlo.rs for the
//! asserted, paper-independent invariants.

#[cfg(not(feature = "validation"))]
fn main() {
    eprintln!("monte_carlo: rebuild with `--features validation` (needs rand, csv, plotters).");
}

#[cfg(feature = "validation")]
fn main() {
    harness::main();
}

#[cfg(feature = "validation")]
mod harness {
    use koenig_planner::cost::Piecewise;
    use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
    use std::f64::consts::TAU;

    /// Chief semimajor axis a_c [m] — the I/O scaling factor (spec §5.5).
    pub const A_C: f64 = 25_000e3;
    /// Per-ROE Gaussian std, metre-scaled (σ = 1 km; spec §6 Phase 6).
    #[allow(dead_code)] // used in later tasks (sweep samplers)
    pub const SIGMA_M: f64 = 1000.0;
    /// Documented constant seed (portable StdRng) — "koenig" in hex-ish.
    pub const SEED: u64 = 0x6F_656E_6967;
    /// Worked-example window [s].
    pub const T_I: f64 = 0.0;
    pub const T_F: f64 = 117_990.0;
    /// Fig. 8 grid step [s] (Table III 30 s grid → 3934 candidate times).
    pub const GRID_DT: f64 = 30.0;

    /// Table III chief mean absolute orbit (angles in radians).
    pub fn worked_example_chief() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            A_C,
            0.7,
            40.0_f64.to_radians(),
            358.0_f64.to_radians(),
            0.0,
            180.0_f64.to_radians(),
        )
    }

    /// J2 mean-ROE dynamics for the worked-example window.
    pub fn worked_example_dynamics() -> J2Roe {
        J2Roe::new(worked_example_chief(), T_I, T_F)
    }

    /// eq. 49 piecewise cost (FaceMax in 2-hr perigee windows, Norm2 elsewhere).
    pub fn worked_example_cost() -> Piecewise {
        Piecewise::new(TAU / worked_example_chief().mean_motion())
    }

    pub fn main() {
        let dynamics = worked_example_dynamics();
        let cost = worked_example_cost();
        let _ = (&dynamics, &cost); // wired into the sweeps in later tasks
        println!("koenig-planner Monte Carlo harness (Phase 6)");
        println!("  problem            : worked example (Table III chief, eq. 49 cost)");
        println!("  window [s]         : [{T_I}, {T_F}]");
        println!(
            "  mean motion [rad/s]: {:.6e}",
            worked_example_chief().mean_motion()
        );
        println!("  Fig. 8 grid        : dt = {GRID_DT} s");
        println!("  seed               : {SEED:#x}");
    }
}
