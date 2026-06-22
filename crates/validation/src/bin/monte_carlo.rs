//! Monte-Carlo validation harness CLI — reproduces [KD20] Fig. 7/8/9.
//! Run: cargo run --release -p koenig-damico-planner-validation --features figures \
//!        --bin monte_carlo [fig7|fig8|fig9]   (default: all)

use koenig_damico_planner_validation::{
    fig7, fig8, fig9, worked_example_cost, worked_example_dynamics, SEED,
};

fn main() {
    let which = std::env::args().nth(1);
    if let Some(arg) = which.as_deref() {
        if !matches!(arg, "fig7" | "fig8" | "fig9") {
            eprintln!("usage: monte_carlo [fig7|fig8|fig9]   (default: all)");
            std::process::exit(2);
        }
    }
    std::fs::create_dir_all("target").ok();
    let dynamics = worked_example_dynamics();
    let cost = worked_example_cost();
    println!("koenig-damico-planner Monte Carlo harness  seed={SEED:#x}");
    let all = which.is_none();
    if all || which.as_deref() == Some("fig7") {
        fig7(&dynamics, &cost);
    }
    if all || which.as_deref() == Some("fig8") {
        fig8(&dynamics, &cost);
    }
    if all || which.as_deref() == Some("fig9") {
        fig9(&dynamics, &cost);
    }
}
