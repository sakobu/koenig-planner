//! Worked example (paper §7 / Table III) — run with the FD-verified J2 mean-ROE
//! dynamics. Reports the planner's computed maneuver plan plus the self-consistent
//! dual lower bound (the exact all-times SOCP optimum), and the Fig. 7 contact
//! curve. Run:
//!
//!   cargo run --example mdot
//!   cargo run --example mdot --features validation   # also writes the Fig. 7 CSV
//!
//! NOTE on reproducing Koenig & D'Amico's published numbers: the paper reports a
//! 3-maneuver, 82.4 mm/s plan. Our dynamics are independently finite-difference
//! verified (tests/fd_stm.rs, tests/fd_b_matrix.rs) at this exact orbit, and in
//! these dynamics the paper's own Table IV maneuvers do NOT reconstruct the
//! Table III target (≈65% residual), while the true optimum is ≈80.9 mm/s. The
//! paper's printed STM also contains a transcription typo we had to correct (the
//! delta-lambda secular drift was -1.5 n dt^2 instead of the dimensionally-correct
//! -1.5 n dt). So this example validates the *math* (FD-verified dynamics + a
//! self-consistent dual/primal), not bit-for-bit reproduction of the paper's
//! (internally inconsistent) worked-example figures.

use koenig_planner::cost::Piecewise;
use koenig_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_planner::{refine_socp, solve, CostModel, Dynamics, Pseudostate, SolveParams, TimeGrid};
use std::f64::consts::TAU;

/// Chief semimajor axis a_c [m] — the I/O scaling factor (spec §5.5).
const A_C: f64 = 25_000e3;
/// Table III target pseudostate in metres (= a_c * w_nd).
const W_METRES: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

fn main() {
    // --- Table III inputs. ---
    let chief = AbsoluteOrbit::new(
        A_C,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let (t_i, t_f) = (0.0, 117_990.0);
    let dynamics = J2Roe::new(chief, t_i, t_f);
    let cost = Piecewise::new(TAU / chief.mean_motion()); // eq.49 perigee windows
    let w = Pseudostate::from_row_slice(&W_METRES) / A_C; // dimensionless w_nd
    let grid = TimeGrid::uniform(t_i, t_f, 30.0); // 3934 candidate times
    let params = SolveParams::default();

    let sol = solve(&dynamics, &cost, w, grid, &params).expect("worked example should solve");

    // --- Exact discretized dual (all-times SOCP) = the self-consistent optimum. ---
    let rows: Vec<_> = grid
        .times()
        .map(|t| cost.at(t).cone_constraints(&dynamics.gamma(t)))
        .collect();
    let exact_dual = refine_socp(&w, &rows).expect("exact SOCP").objective;

    // --- Report. ---
    println!("Koenig planner — worked example (Table III), FD-verified dynamics");
    println!("  candidate times      : {}", grid.len());
    println!("  iterations           : {}", sol.iterations);
    println!("  residual w_err/w     : {:.3e}", sol.residual);
    println!("  total dv (computed)  : {:.4} mm/s", sol.total_dv * 1e3);
    println!(
        "  dual lower bound     : {:.4} mm/s  (exact all-times SOCP optimum)",
        exact_dual * 1e3
    );
    println!(
        "  refinement objective : {:.4} mm/s",
        sol.lambda.dot(&w) * 1e3
    );
    println!("  paper (Koenig 7)     : 82.4 mm/s primal / 82.0 mm/s bound");
    println!("  maneuvers ({}):", sol.maneuvers.len());
    println!(
        "    {:>10}  {:>9}  {:>9}  {:>9}  {:>9}",
        "t [s]", "u_R", "u_T", "u_N", "|dv|"
    );
    for m in &sol.maneuvers {
        println!(
            "    {:>10.0}  {:>9.3}  {:>9.3}  {:>9.3}  {:>9.4}",
            m.t,
            m.dv[0] * 1e3,
            m.dv[1] * 1e3,
            m.dv[2] * 1e3,
            m.dv.norm() * 1e3
        );
    }

    // --- Fig. 7: contact function g(t) = g_{U(1,t)}(Gamma^T(t) lambda) over the grid. ---
    let curve: Vec<(f64, f64)> = grid
        .times()
        .map(|t| {
            let y = dynamics.gamma(t).transpose() * sol.lambda;
            (t, cost.at(t).contact(y))
        })
        .collect();
    let max_g = curve
        .iter()
        .map(|&(_, g)| g)
        .fold(f64::NEG_INFINITY, f64::max);
    println!(
        "  Fig. 7 max_t g       : {:.6}  (<= 1 + eps_cost = {:.2})",
        max_g,
        1.0 + params.eps_cost
    );

    #[cfg(feature = "validation")]
    {
        let path = "target/fig7_contact.csv";
        let mut wtr = csv::Writer::from_path(path).expect("open fig7 csv");
        wtr.write_record(["t_s", "g"]).expect("header");
        for (t, g) in &curve {
            wtr.write_record(&[t.to_string(), g.to_string()])
                .expect("row");
        }
        wtr.flush().expect("flush");
        println!(
            "  Fig. 7 curve         : written to {path} ({} rows)",
            curve.len()
        );
    }
    #[cfg(not(feature = "validation"))]
    println!("  Fig. 7 curve         : (build with --features validation to write target/fig7_contact.csv)");

    // --- Self-checks (what is actually true of the FD-verified pipeline). ---
    assert!(
        sol.iterations >= 1 && sol.iterations <= 50,
        "did not converge"
    );
    assert!(
        max_g <= 1.0 + params.eps_cost + 1e-6,
        "max_t g exceeds tolerance"
    );
    // The refinement finds the true discretized dual optimum (self-consistency).
    let rel = (sol.lambda.dot(&w) - exact_dual).abs() / exact_dual;
    assert!(
        rel < 1e-2,
        "refinement dual {} vs exact {}",
        sol.lambda.dot(&w),
        exact_dual
    );
}
