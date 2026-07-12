//! Criterion benches for the reachable-set sweep hot path.
//!
//!   cargo bench --bench sweep_dual
//!
//! Two groups, all on the public API (no clarabel-internal coupling):
//!
//! * `point-cost-HEO` — per-target cost on the KD20 Table III worked example:
//!   one `refine_socp` (exact dual) and one full `solve()` pipeline (Alg 1→2→3),
//!   at a coarse (~256-pt) and the full (3934-pt) grid, plus the Γ-row assembly
//!   alone. The `refine_socp` vs `assembly` pair is the assembly-vs-solve split;
//!   the `solve/*` pair is the point-click tier.
//!
//! * `trace-256dir` — one 256-direction reachable-set boundary trace at the coarse
//!   landscape grid, `sweep_dual` vs looping `solve()`, across all four demo presets
//!   (`crates/wasm/www/src/lib/defaults.ts`). This is the landscape-engine decision.
//!
//! Not here (they are validations / one-shot diagnostics, not throughput benches):
//! the within-clarabel setup-vs-solve split (needs a clarabel replica), the
//! `solve` vs `sweep_dual` c* agreement, and the Q2 `SolveQuality` / Q3 burn-count
//! demonstrations. Those belong in `tests/` once the signals they exercise land.

use std::f64::consts::TAU;
use std::hint::black_box;
use std::time::Duration;

use criterion::measurement::WallTime;
use criterion::{
    criterion_group, criterion_main, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};

use koenig_damico_planner::cost::{Norm2, Piecewise};
use koenig_damico_planner::dynamics::{AbsoluteOrbit, J2Roe};
use koenig_damico_planner::solver::{refine_socp, sweep_dual};
use koenig_damico_planner::types::ConicRows;
use koenig_damico_planner::{
    solve, CostModel, Dynamics, Pseudostate, SolveParams, SublevelSet, TimeGrid,
};

const A_C: f64 = 25_000e3;
/// KD20 Table III target pseudostate `[δa, δλ, δe_x, δe_y, δi_x, δi_y]` (meters).
const W_HEO: [f64; 6] = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0];

/// Constant Norm2 cost (the presets' `norm2`; mirrors api's internal `ConstNorm2`).
struct ConstNorm2(Norm2);
impl CostModel for ConstNorm2 {
    fn at(&self, _t: f64) -> &dyn SublevelSet {
        &self.0
    }
}

/// Assemble the per-time conic rows (Γ(t) propagation + `cone_constraints`), exactly
/// as `sweep_dual` does once per window.
fn build_rows<C: CostModel>(dyn_: &J2Roe, cost: &C, grid: &TimeGrid) -> Vec<ConicRows> {
    grid.times()
        .map(|t| cost.at(t).cone_constraints(&dyn_.gamma(t).unwrap()))
        .collect()
}

/// 256 unit directions in the (δe_x, δi_y) plane — mixes an in-plane and a
/// cross-track axis, so it also exercises GEO's near-equatorial B(t) conditioning.
fn directions(n: usize) -> Vec<Pseudostate> {
    (0..n)
        .map(|k| {
            let th = TAU * (k as f64) / (n as f64);
            let mut d = Pseudostate::zeros();
            d[2] = th.cos();
            d[5] = th.sin();
            d
        })
        .collect()
}

fn heo_chief() -> AbsoluteOrbit {
    AbsoluteOrbit::new(
        A_C,
        0.7,
        40f64.to_radians(),
        358f64.to_radians(),
        0.0,
        180f64.to_radians(),
    )
}

/// Add the `sweep_dual` vs `loop solve()` pair for one preset to the trace group.
fn add_trace<C: CostModel>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    preset: &str,
    dyn_: &J2Roe,
    cost: &C,
    coarse: &TimeGrid,
) {
    let dirs = directions(256);
    let sp = SolveParams::default();
    group.throughput(Throughput::Elements(dirs.len() as u64));
    group.bench_function(BenchmarkId::new("sweep_dual", preset), |b| {
        b.iter(|| {
            let r = sweep_dual(dyn_, cost, coarse, black_box(&dirs)).unwrap();
            black_box(
                r.iter()
                    .filter(|p| p.feasible)
                    .map(|p| p.c_star)
                    .sum::<f64>(),
            )
        })
    });
    group.bench_function(BenchmarkId::new("loop_solve", preset), |b| {
        b.iter(|| {
            black_box(
                dirs.iter()
                    .filter_map(|d| solve(dyn_, cost, *d, *coarse, &sp).ok())
                    .map(|s| s.total_dv)
                    .sum::<f64>(),
            )
        })
    });
}

fn bench_point_cost(c: &mut Criterion) {
    let chief = heo_chief();
    let dyn_ = J2Roe::new(chief, 0.0, 117_990.0).unwrap();
    let cost = Piecewise::new(TAU / chief.mean_motion()).unwrap();
    let w = Pseudostate::from_row_slice(&W_HEO) / A_C;
    let sp = SolveParams::default();

    let coarse = TimeGrid::uniform(0.0, 117_990.0, 461.0).unwrap(); // ~256 pts
    let full = TimeGrid::uniform(0.0, 117_990.0, 30.0).unwrap(); // 3934 pts
    let rows_coarse = build_rows(&dyn_, &cost, &coarse);
    let rows_full = build_rows(&dyn_, &cost, &full);

    let mut g = c.benchmark_group("point-cost-HEO");
    g.sample_size(30); // one full refine_socp is ~80 ms; keep the group bounded.

    g.bench_function("refine_socp/coarse-256", |b| {
        b.iter(|| {
            refine_socp(black_box(&w), black_box(&rows_coarse))
                .unwrap()
                .objective
        })
    });
    g.bench_function("refine_socp/full-3934", |b| {
        b.iter(|| {
            refine_socp(black_box(&w), black_box(&rows_full))
                .unwrap()
                .objective
        })
    });
    g.bench_function("solve/coarse-256", |b| {
        b.iter(|| solve(&dyn_, &cost, w, coarse, &sp).unwrap().total_dv)
    });
    g.bench_function("solve/full-3934", |b| {
        b.iter(|| solve(&dyn_, &cost, w, full, &sp).unwrap().total_dv)
    });
    // Γ-row assembly alone: the "assembly" half of the assembly-vs-solve split.
    // Fence the whole `rows` value, not just its length, so the `cone_constraints`
    // work (the bulk of the assembly cost) can't be elided as dead code.
    g.bench_function("assembly/coarse-256", |b| {
        b.iter(|| {
            let rows = build_rows(&dyn_, &cost, &coarse);
            black_box(&rows);
            rows.len()
        })
    });
    g.finish();
}

fn bench_trace(c: &mut Criterion) {
    let mut g = c.benchmark_group("trace-256dir");
    // One trace is ~0.1–1.2 s; shorten warm-up and cap samples so `cargo bench`
    // stays reasonable while still statistically meaningful.
    g.sample_size(10);
    g.warm_up_time(Duration::from_secs(1));

    // 1 — HEO worked example (piecewise). Native 3934 pts is too slow for a 256-dir
    // sweep (~20 s), so the trace runs on the coarse landscape grid.
    let heo = heo_chief();
    let heo_dyn = J2Roe::new(heo, 0.0, 117_990.0).unwrap();
    let heo_cost = Piecewise::new(TAU / heo.mean_motion()).unwrap();
    add_trace(
        &mut g,
        "HEO-piecewise",
        &heo_dyn,
        &heo_cost,
        &TimeGrid::uniform(0.0, 117_990.0, 461.0).unwrap(),
    );

    // 2 — LEO sun-sync formation (norm2).
    let sso = AbsoluteOrbit::new(
        6_900e3,
        0.001,
        97.485f64.to_radians(),
        90f64.to_radians(),
        0.0,
        0.0,
    );
    let sso_dyn = J2Roe::new(sso, 0.0, 11_400.0).unwrap();
    add_trace(
        &mut g,
        "LEO-sso-norm2",
        &sso_dyn,
        &ConstNorm2(Norm2),
        &TimeGrid::uniform(0.0, 11_400.0, 45.0).unwrap(),
    );

    // 3 — LEO co-elliptic hold (norm2).
    let coe = AbsoluteOrbit::new(6_878e3, 0.0005, 51.6f64.to_radians(), 0.0, 0.0, 0.0);
    let coe_dyn = J2Roe::new(coe, 0.0, 11_400.0).unwrap();
    add_trace(
        &mut g,
        "LEO-coelliptic-norm2",
        &coe_dyn,
        &ConstNorm2(Norm2),
        &TimeGrid::uniform(0.0, 11_400.0, 45.0).unwrap(),
    );

    // 4 — GEO relative slot offset (norm2). Near-equatorial (i=0.05°): B(t)
    // cross-track terms are ill-conditioned — the convergence stress case.
    let geo = AbsoluteOrbit::new(42_164e3, 0.0002, 0.05f64.to_radians(), 0.0, 0.0, 0.0);
    let geo_dyn = J2Roe::new(geo, 0.0, 43_082.0).unwrap();
    add_trace(
        &mut g,
        "GEO-slot-norm2",
        &geo_dyn,
        &ConstNorm2(Norm2),
        &TimeGrid::uniform(0.0, 43_082.0, 170.0).unwrap(),
    );

    g.finish();
}

criterion_group!(benches, bench_point_cost, bench_trace);
criterion_main!(benches);
