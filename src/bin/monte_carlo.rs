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
//! Fig. 8 reproduces the paper's THREE initialization schemes (Koenig & D'Amico 2020,
//! p.11): n=2 seeds only the window endpoints {t_i, t_f} (deliberate worst case), n=6
//! is Algorithm 1 (the six largest-contact times), n=10 is ten evenly-spaced times —
//! NOT a single `n_init` count knob. Validation stance (spec §6 Phase 6): we REPORT our
//! iteration distributions and compare to the paper's 4.90/3.99/3.31 means as
//! *reference*, not as a pass/fail target. See tests/monte_carlo.rs for the asserted,
//! paper-independent invariants.

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
    use koenig_planner::{
        solve, solve_from_initial_times, CostModel, Dynamics, PlannerError, Pseudostate, Solution,
        SolveParams, TimeGrid,
    };
    use plotters::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal};
    use std::f64::consts::TAU;
    use std::time::Instant;

    /// Fig. 8 sample count (paper: 200).
    pub const N_MC: usize = 200;
    /// The three Fig. 8 initialization schemes (Koenig & D'Amico 2020, p.11),
    /// labelled by candidate count. These are three DISTINCT seedings, not a single
    /// `n_init` count knob — see [`InitScheme`]. Routing all three through Algorithm 1
    /// would make the n=2/n=10 columns measure a different (stronger) seed than the
    /// paper's, so each is dispatched to its own seeding.
    pub const FIG8_SCHEMES: [(usize, InitScheme); 3] = [
        (2, InitScheme::Endpoints),
        (6, InitScheme::LargestG(6)),
        (10, InitScheme::EvenlySpaced(10)),
    ];
    /// Paper's reported mean iterations for the three schemes (reference, not a target).
    pub const PAPER_MEANS: [f64; 3] = [4.90, 3.99, 3.31];

    /// Chief semimajor axis a_c `[m]` — the I/O scaling factor (spec §5.5).
    pub const A_C: f64 = 25_000e3;
    /// Per-ROE Gaussian std, metre-scaled (σ = 1 km; spec §6 Phase 6).
    pub const SIGMA_M: f64 = 1000.0;
    /// Documented constant seed (portable StdRng) — "koenig" in hex-ish.
    pub const SEED: u64 = 0x6F_656E_6967;
    /// Worked-example window `[s]`.
    pub const T_I: f64 = 0.0;
    pub const T_F: f64 = 117_990.0;
    /// Fig. 8 grid step `[s]` (Table III 30 s grid → 3934 candidate times).
    pub const GRID_DT: f64 = 30.0;
    /// Fig. 9 grid sizes (10 → 10⁶). 10⁶ is ~150 MB Γ cache / multi-second; documented.
    pub const FIG9_SIZES: [usize; 6] = [10, 100, 1_000, 10_000, 100_000, 1_000_000];

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
            .expect("worked-example chief is a valid elliptic, inclined orbit")
    }

    /// eq. 49 piecewise cost (FaceMax in 2-hr perigee windows, Norm2 elsewhere).
    pub fn worked_example_cost() -> Piecewise {
        Piecewise::new(TAU / worked_example_chief().mean_motion())
    }

    pub fn main() {
        let which = std::env::args().nth(1);
        if let Some(arg) = which.as_deref() {
            if arg != "fig8" && arg != "fig9" {
                eprintln!("usage: monte_carlo [fig8|fig9]   (default: both)");
                std::process::exit(2);
            }
        }
        std::fs::create_dir_all("target").ok();
        let dynamics = worked_example_dynamics();
        let cost = worked_example_cost();

        println!("koenig-planner Monte Carlo harness (Phase 6)  seed={SEED:#x}");
        let run_8 = matches!(which.as_deref(), None | Some("fig8"));
        if run_8 {
            fig8(&dynamics, &cost);
        }
        let run_9 = matches!(which.as_deref(), None | Some("fig9"));
        if run_9 {
            fig9(&dynamics, &cost);
        }
    }

    /// `n` random target pseudostates as dimensionless `w_nd`: each of the 6 ROE
    /// components `~ Normal(0, σ = SIGMA_M metres)`, then divided by `a_c`
    /// (spec §6 Phase 6 sampling convention). `StdRng` is portable, so a fixed
    /// `seed` yields identical samples on every platform.
    pub fn sample_pseudostates(n: usize, seed: u64) -> Vec<Pseudostate> {
        let mut rng = StdRng::seed_from_u64(seed);
        let normal = Normal::new(0.0_f64, SIGMA_M).expect("σ > 0 is a valid normal");
        (0..n)
            .map(|_| {
                let mut comp = [0.0_f64; 6];
                for c in comp.iter_mut() {
                    *c = normal.sample(&mut rng);
                }
                Pseudostate::from_row_slice(&comp) / A_C
            })
            .collect()
    }

    /// A Fig. 8 initialization scheme — the paper used three DISTINCT seedings for
    /// the 2/6/10-candidate columns (Koenig & D'Amico 2020, "Sensitivity to Poor
    /// Initialization", p.11), not one `n_init` knob.
    #[derive(Clone, Copy)]
    pub enum InitScheme {
        /// Algorithm 1: the `n` largest-contact coarse times (the nominal 6-time seed).
        LargestG(usize),
        /// Only the window endpoints `{t_i, t_f}` — the paper's worst-case 2-time seed.
        Endpoints,
        /// `n` times evenly spaced over `[t_i, t_f]` — the paper's 10-time seed.
        EvenlySpaced(usize),
    }

    /// Short label for the CSV / printout.
    fn scheme_name(scheme: InitScheme) -> &'static str {
        match scheme {
            InitScheme::LargestG(_) => "largest_g",
            InitScheme::Endpoints => "endpoints",
            InitScheme::EvenlySpaced(_) => "evenly_spaced",
        }
    }

    /// `n` times evenly spaced over `[t_i, t_f]`, inclusive of both endpoints.
    fn evenly_spaced_times(t_i: f64, t_f: f64, n: usize) -> Vec<f64> {
        if n <= 1 {
            return vec![t_i];
        }
        (0..n)
            .map(|k| t_i + (k as f64) * (t_f - t_i) / ((n - 1) as f64))
            .collect()
    }

    /// Solve one Fig. 8 sample under a given initialization scheme on the fixed 30 s
    /// grid. `LargestG` uses Algorithm 1 (`solve`); the other two seed refinement
    /// directly with explicit times (`solve_from_initial_times`), as the paper does.
    fn solve_scheme<D: Dynamics, C: CostModel>(
        dynamics: &D,
        cost: &C,
        w: Pseudostate,
        grid: TimeGrid,
        scheme: InitScheme,
    ) -> Result<Solution, PlannerError> {
        match scheme {
            InitScheme::LargestG(n) => {
                let params = SolveParams {
                    n_init: n,
                    ..SolveParams::default()
                };
                solve(dynamics, cost, w, grid, &params)
            }
            InitScheme::Endpoints => solve_from_initial_times(
                dynamics,
                cost,
                w,
                grid,
                &SolveParams::default(),
                &[grid.t_i, grid.t_f],
            ),
            InitScheme::EvenlySpaced(n) => {
                let times = evenly_spaced_times(grid.t_i, grid.t_f, n);
                solve_from_initial_times(dynamics, cost, w, grid, &SolveParams::default(), &times)
            }
        }
    }

    /// One Fig. 8 sample outcome.
    #[derive(Clone, Copy)]
    pub struct Fig8Row {
        pub n_init: usize,
        pub scheme: &'static str,
        pub sample: usize,
        pub iterations: usize,
        pub residual: f64,
        pub total_dv: f64,
    }

    /// Per-`n_init` summary statistics.
    pub struct Fig8Stat {
        pub n_init: usize,
        pub n: usize,
        pub mean_iters: f64,
        pub frac_le8: f64,
        pub max_iters: usize,
        pub max_residual: f64,
    }

    /// Run every `(scheme, w)` pair on the fixed 30 s grid; collect per-sample
    /// outcomes and a count of solver failures (expect 0). Each scheme uses the
    /// paper's own seeding (endpoints / largest-g / evenly-spaced), so the columns
    /// are comparable to Fig. 8.
    pub fn run_fig8<D: Dynamics, C: CostModel>(
        dynamics: &D,
        cost: &C,
        ws: &[Pseudostate],
        schemes: &[(usize, InitScheme)],
    ) -> (Vec<Fig8Row>, usize) {
        let grid = TimeGrid::uniform(T_I, T_F, GRID_DT).expect("valid sweep grid");
        let mut rows = Vec::with_capacity(ws.len() * schemes.len());
        let mut failures = 0usize;
        for &(n_init, scheme) in schemes {
            let name = scheme_name(scheme);
            for (sample, &w) in ws.iter().enumerate() {
                match solve_scheme(dynamics, cost, w, grid, scheme) {
                    Ok(sol) => rows.push(Fig8Row {
                        n_init,
                        scheme: name,
                        sample,
                        iterations: sol.iterations,
                        residual: sol.residual,
                        total_dv: sol.total_dv,
                    }),
                    Err(_) => failures += 1,
                }
            }
        }
        (rows, failures)
    }

    /// Group rows by scheme label (not the `n_init` count, which could collide
    /// across future schemes) and compute mean iterations, fraction ≤ 8
    /// iterations, max iterations, and max residual.
    pub fn summarize_fig8(rows: &[Fig8Row], schemes: &[(usize, InitScheme)]) -> Vec<Fig8Stat> {
        debug_assert!(
            {
                let mut names: Vec<&str> = schemes.iter().map(|&(_, s)| scheme_name(s)).collect();
                names.sort_unstable();
                names.dedup();
                names.len() == schemes.len()
            },
            "Fig.8 scheme labels must be unique to group correctly"
        );
        schemes
            .iter()
            .map(|&(n_init, scheme)| {
                let name = scheme_name(scheme);
                let group: Vec<&Fig8Row> = rows.iter().filter(|r| r.scheme == name).collect();
                let n = group.len();
                let denom = n.max(1) as f64;
                let mean_iters = group.iter().map(|r| r.iterations as f64).sum::<f64>() / denom;
                let frac_le8 = group.iter().filter(|r| r.iterations <= 8).count() as f64 / denom;
                let max_iters = group.iter().map(|r| r.iterations).max().unwrap_or(0);
                let max_residual = group.iter().map(|r| r.residual).fold(0.0, f64::max);
                Fig8Stat {
                    n_init,
                    n,
                    mean_iters,
                    frac_le8,
                    max_iters,
                    max_residual,
                }
            })
            .collect()
    }

    /// Write the per-sample Fig. 8 rows to `path` as CSV.
    pub fn write_fig8_csv(path: &str, rows: &[Fig8Row]) -> csv::Result<()> {
        let mut w = csv::Writer::from_path(path)?;
        w.write_record([
            "n_init",
            "scheme",
            "sample",
            "iterations",
            "residual",
            "total_dv",
        ])?;
        for r in rows {
            w.write_record(&[
                r.n_init.to_string(),
                r.scheme.to_string(),
                r.sample.to_string(),
                r.iterations.to_string(),
                format!("{:.6e}", r.residual),
                format!("{:.9e}", r.total_dv),
            ])?;
        }
        w.flush()?;
        Ok(())
    }

    /// Fig. 8 driver: sample, sweep, summarize, print, and write the CSV.
    pub fn fig8<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C) {
        let ws = sample_pseudostates(N_MC, SEED);
        let (rows, failures) = run_fig8(dynamics, cost, &ws, &FIG8_SCHEMES);
        let stats = summarize_fig8(&rows, &FIG8_SCHEMES);

        println!("\nFig. 8 — Algorithm-2 iteration distribution ({N_MC} samples/scheme)");
        println!("  paper's three seedings: n=2 endpoints, n=6 largest-g, n=10 evenly-spaced");
        println!(
            "  {:>6}  {:>14}  {:>5}  {:>10}  {:>8}  {:>11}  {:>12}",
            "n_init", "scheme", "n", "mean_iters", "frac<=8", "max_iters", "max_residual"
        );
        for ((s, &(_, scheme)), paper) in stats
            .iter()
            .zip(FIG8_SCHEMES.iter())
            .zip(PAPER_MEANS.iter())
        {
            println!(
                "  {:>6}  {:>14}  {:>5}  {:>10.3}  {:>8.3}  {:>11}  {:>12.2e}   (paper {:.2})",
                s.n_init,
                scheme_name(scheme),
                s.n,
                s.mean_iters,
                s.frac_le8,
                s.max_iters,
                s.max_residual,
                paper
            );
        }
        if failures > 0 {
            println!("  WARNING: {failures} solve(s) failed (expected 0).");
        }

        let path = "target/fig8_iterations.csv";
        match write_fig8_csv(path, &rows) {
            Ok(()) => println!("  rows written         : {path} ({} rows)", rows.len()),
            Err(e) => eprintln!("  CSV write failed     : {e}"),
        }

        let cdf: Vec<(usize, Vec<(f64, f64)>)> = FIG8_SCHEMES
            .iter()
            .map(|&(n_init, scheme)| {
                let name = scheme_name(scheme);
                let counts: Vec<usize> = rows
                    .iter()
                    .filter(|r| r.scheme == name)
                    .map(|r| r.iterations)
                    .collect();
                (n_init, empirical_cdf(&counts))
            })
            .collect();
        match plot_fig8_cdf("target/fig8_cdf.png", &cdf) {
            Ok(()) => println!("  CDF plot             : target/fig8_cdf.png"),
            Err(e) => eprintln!("  PNG write failed     : {e}"),
        }
    }

    /// One Fig. 9 timing outcome.
    #[derive(Clone, Copy)]
    pub struct Fig9Row {
        pub grid_len: usize,
        pub dt: f64,
        pub seconds: f64,
        pub iterations: usize,
        pub residual: f64,
    }

    /// Time `solve` once per requested grid size on the fixed window. For each size:
    /// `dt = (t_f - t_i)/(n-1)`, one warmup solve (discarded), one timed solve. The
    /// actual `grid.len()` is recorded (it may differ from `n` by ±1 due to rounding
    /// in `TimeGrid::len`). Timing shape is `w`-independent; use a single fixed `w`.
    pub fn run_fig9<D: Dynamics, C: CostModel>(
        dynamics: &D,
        cost: &C,
        w: Pseudostate,
        sizes: &[usize],
    ) -> Vec<Fig9Row> {
        let params = SolveParams::default();
        sizes
            .iter()
            .map(|&n| {
                let dt = (T_F - T_I) / (n.max(2) - 1) as f64;
                let grid = TimeGrid::uniform(T_I, T_F, dt).expect("valid sweep grid");
                let _ = solve(dynamics, cost, w, grid, &params); // warmup
                let start = Instant::now();
                let result = solve(dynamics, cost, w, grid, &params);
                let seconds = start.elapsed().as_secs_f64();
                match result {
                    Ok(s) => Fig9Row {
                        grid_len: grid.len(),
                        dt,
                        seconds,
                        iterations: s.iterations,
                        residual: s.residual,
                    },
                    Err(_) => Fig9Row {
                        grid_len: grid.len(),
                        dt,
                        seconds,
                        iterations: 0,
                        residual: f64::NAN,
                    },
                }
            })
            .collect()
    }

    /// Write the Fig. 9 timing rows to `path` as CSV.
    pub fn write_fig9_csv(path: &str, rows: &[Fig9Row]) -> csv::Result<()> {
        let mut w = csv::Writer::from_path(path)?;
        w.write_record(["grid_len", "dt_s", "seconds", "iterations", "residual"])?;
        for r in rows {
            w.write_record(&[
                r.grid_len.to_string(),
                format!("{:.6e}", r.dt),
                format!("{:.6e}", r.seconds),
                r.iterations.to_string(),
                format!("{:.3e}", r.residual),
            ])?;
        }
        w.flush()?;
        Ok(())
    }

    /// Fig. 9 driver: time the Table III `w` across `FIG9_SIZES`, print, write CSV.
    pub fn fig9<D: Dynamics, C: CostModel>(dynamics: &D, cost: &C) {
        let w = sample_pseudostates(1, SEED)[0];
        println!("\nFig. 9 — solve time vs |T| (10⁶ is multi-second / ~150 MB)");
        let rows = run_fig9(dynamics, cost, w, &FIG9_SIZES);
        println!(
            "  {:>10}  {:>12}  {:>10}  {:>6}  {:>10}",
            "grid_len", "dt_s", "seconds", "iters", "residual"
        );
        for r in &rows {
            println!(
                "  {:>10}  {:>12.4e}  {:>10.4}  {:>6}  {:>10.2e}",
                r.grid_len, r.dt, r.seconds, r.iterations, r.residual
            );
        }
        let path = "target/fig9_timing.csv";
        match write_fig9_csv(path, &rows) {
            Ok(()) => println!("  rows written         : {path} ({} rows)", rows.len()),
            Err(e) => eprintln!("  CSV write failed     : {e}"),
        }

        match plot_fig9_timing("target/fig9_timing.png", &rows) {
            Ok(()) => println!("  timing plot          : target/fig9_timing.png"),
            Err(e) => eprintln!("  PNG write failed     : {e}"),
        }
    }

    /// Empirical CDF of iteration counts as `(value, fraction ≤ value)` over the
    /// distinct sorted values, anchored at `(min-1, 0)` so the curve starts at 0.
    pub fn empirical_cdf(counts: &[usize]) -> Vec<(f64, f64)> {
        assert!(
            !counts.is_empty(),
            "empirical_cdf needs at least one sample"
        );
        let n = counts.len() as f64;
        let mut sorted = counts.to_vec();
        sorted.sort_unstable();
        let mut pts = vec![((*sorted.first().unwrap() as f64) - 1.0, 0.0)];
        let mut i = 0;
        while i < sorted.len() {
            let v = sorted[i];
            let mut j = i;
            while j < sorted.len() && sorted[j] == v {
                j += 1;
            }
            pts.push((v as f64, j as f64 / n));
            i = j;
        }
        pts
    }

    /// Plot one step-CDF per `n_init` to a PNG.
    pub fn plot_fig8_cdf(
        path: &str,
        series: &[(usize, Vec<(f64, f64)>)],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (900, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let x_max = series
            .iter()
            .flat_map(|(_, p)| p.iter().map(|&(x, _)| x))
            .fold(1.0_f64, f64::max);
        let mut chart = ChartBuilder::on(&root)
            .caption("Fig. 8 - Algorithm-2 iteration CDF", ("sans-serif", 28))
            .margin(12)
            .x_label_area_size(45)
            .y_label_area_size(55)
            .build_cartesian_2d(0f64..(x_max + 1.0), 0f64..1.02f64)?;
        chart
            .configure_mesh()
            .x_desc("iterations")
            .y_desc("empirical CDF")
            .draw()?;
        let palette = [RED, BLUE, GREEN];
        for (i, (n_init, pts)) in series.iter().enumerate() {
            let color = palette[i % palette.len()];
            chart
                .draw_series(LineSeries::new(pts.iter().cloned(), color.stroke_width(2)))?
                .label(format!("n_init = {n_init}"))
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], color));
        }
        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.85))
            .border_style(BLACK)
            .draw()?;
        root.present()?;
        Ok(())
    }

    /// Plot solve-time vs |T| on log-log axes to a PNG.
    pub fn plot_fig9_timing(
        path: &str,
        rows: &[Fig9Row],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (900, 600)).into_drawing_area();
        root.fill(&WHITE)?;
        let xs: Vec<f64> = rows.iter().map(|r| r.grid_len as f64).collect();
        let ys: Vec<f64> = rows.iter().map(|r| r.seconds.max(1e-6)).collect();
        let x_lo = xs.iter().cloned().fold(f64::INFINITY, f64::min).max(1.0);
        let x_hi = xs.iter().cloned().fold(0.0, f64::max).max(10.0);
        let y_lo = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_hi = ys.iter().cloned().fold(0.0, f64::max);
        let mut chart = ChartBuilder::on(&root)
            .caption("Fig. 9 - solve time vs |T|", ("sans-serif", 28))
            .margin(12)
            .x_label_area_size(45)
            .y_label_area_size(70)
            .build_cartesian_2d(
                (x_lo..x_hi * 1.5).log_scale(),
                (y_lo * 0.5..y_hi * 2.0).log_scale(),
            )?;
        chart
            .configure_mesh()
            .x_desc("|T| (grid size)")
            .y_desc("solve time [s]")
            .draw()?;
        chart.draw_series(LineSeries::new(
            xs.iter().cloned().zip(ys.iter().cloned()),
            BLUE.stroke_width(2),
        ))?;
        chart.draw_series(
            xs.iter()
                .cloned()
                .zip(ys.iter().cloned())
                .map(|(x, y)| Circle::new((x, y), 3, BLUE.filled())),
        )?;
        root.present()?;
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use approx::assert_abs_diff_eq;

        #[test]
        fn sampler_is_deterministic_and_well_scaled() {
            let a = sample_pseudostates(200, SEED);
            let b = sample_pseudostates(200, SEED);
            assert_eq!(a.len(), 200);
            // Determinism: same seed -> identical samples.
            for (x, y) in a.iter().zip(&b) {
                assert_eq!(x, y);
            }
            // Convention: components ~ Normal(0, SIGMA_M / A_C); never near-zero norm.
            let expected_sd = SIGMA_M / A_C;
            let flat: Vec<f64> = a.iter().flat_map(|w| w.iter().copied()).collect();
            let mean = flat.iter().sum::<f64>() / flat.len() as f64;
            let var = flat.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / flat.len() as f64;
            assert_abs_diff_eq!(mean, 0.0, epsilon = expected_sd * 0.15);
            assert_abs_diff_eq!(var.sqrt(), expected_sd, epsilon = expected_sd * 0.15);
            assert!(a.iter().all(|w| w.norm() > 0.0));
        }

        #[test]
        fn different_seeds_differ() {
            let a = sample_pseudostates(8, SEED);
            let b = sample_pseudostates(8, SEED + 1);
            assert!(a.iter().zip(&b).any(|(x, y)| x != y));
        }

        #[test]
        fn fig9_sweep_times_each_size() {
            let dynamics = worked_example_dynamics();
            let cost = worked_example_cost();
            let w = sample_pseudostates(1, SEED)[0];
            let sizes = [10usize, 100];
            let rows = run_fig9(&dynamics, &cost, w, &sizes);
            assert_eq!(rows.len(), 2);
            assert!(rows.iter().all(|r| r.seconds >= 0.0 && r.grid_len >= 2));
            // Finer grid is at least as large in point count.
            assert!(rows[1].grid_len >= rows[0].grid_len);
        }

        #[test]
        fn empirical_cdf_is_monotone_and_ends_at_one() {
            let pts = empirical_cdf(&[3, 3, 4, 5]);
            // Left anchor at 0, then steps at the distinct values, ending at 1.0.
            assert_eq!(pts.first().unwrap().1, 0.0);
            assert!((pts.last().unwrap().1 - 1.0).abs() < 1e-12);
            // (value, fraction <= value): (2,0), (3,0.5), (4,0.75), (5,1.0).
            assert_eq!(pts, vec![(2.0, 0.0), (3.0, 0.5), (4.0, 0.75), (5.0, 1.0)]);
            // Monotone non-decreasing in the fraction.
            assert!(pts.windows(2).all(|w| w[1].1 >= w[0].1));
        }

        #[test]
        fn fig8_sweep_produces_paired_rows_and_low_residual() {
            let dynamics = worked_example_dynamics();
            let cost = worked_example_cost();
            let ws = sample_pseudostates(3, SEED);
            // One Algorithm-1 scheme and the explicit-seed endpoints scheme, so both
            // solve paths (solve / solve_from_initial_times) are exercised.
            let schemes = [
                (2usize, InitScheme::Endpoints),
                (6usize, InitScheme::LargestG(6)),
            ];
            let (rows, failures) = run_fig8(&dynamics, &cost, &ws, &schemes);
            assert_eq!(
                failures, 0,
                "no solve should fail on the worked-example problem"
            );
            assert_eq!(rows.len(), 3 * 2, "one row per (scheme, sample)");
            for r in &rows {
                assert!(
                    r.residual < 1e-3,
                    "row residual {:.3e} too high",
                    r.residual
                );
                assert!((1..=50).contains(&r.iterations));
            }
            let stats = summarize_fig8(&rows, &schemes);
            assert_eq!(stats.len(), 2);
            assert!(stats.iter().all(|s| s.n == 3));
        }
    }
}
