//! Public-API smoke tests for the Phase 0 scaffold.

use approx::assert_abs_diff_eq;
use koenig_planner::{Maneuver, PlannerError, SolveParams, TimeGrid, M, N};
use nalgebra::SVector;

#[test]
fn reexports_are_reachable() {
    assert_eq!(N, 6);
    assert_eq!(M, 3);
}

#[test]
fn maneuver_constructs_and_exposes_fields() {
    let m = Maneuver {
        t: 16050.0,
        dv: SVector::<f64, 3>::new(9.68e-3, -23.02e-3, -25.56e-3),
    };
    assert_abs_diff_eq!(m.t, 16050.0, epsilon = 1e-9);
    assert_eq!(m.dv.len(), 3);
}

#[test]
fn default_params_are_table_iii() {
    let p = SolveParams::default();
    assert_eq!(p.n_init, 6);
    assert_eq!(p.n_coarse, 20);
}

#[test]
fn error_displays_message() {
    let e = PlannerError::InvalidInput("bad w".into());
    assert!(e.to_string().contains("bad w"));
}

#[test]
fn worked_and_hunter_grid_counts() {
    assert_eq!(TimeGrid::uniform(0.0, 117990.0, 30.0).len(), 3934);
    assert_eq!(TimeGrid::uniform(0.0, 39000.0, 10.0).len(), 3901);
}
