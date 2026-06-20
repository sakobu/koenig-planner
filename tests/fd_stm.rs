//! Independent finite-difference verification of the J2 mean-ROE STM Phi(t,t_f).
//!
//! This does NOT use the analytic STM formula. It reconstructs Phi from the eq.51
//! ROE map `x = roe(chief, deputy)` (algebraic) and the eq.50 secular propagation
//! `oe(t_f) = propagate(oe)` (AbsoluteOrbit::propagate), via
//! `Phi = d x(t_f)/d x(t) = J_f * P * J_t^{-1}`, where J_t, J_f are the eq.51
//! Jacobians and P = d propagate / d oe is the secular-propagation Jacobian. Both
//! factors are finite-differenced, so the whole Phi is independent of stm.rs.

use koenig_planner::dynamics::stm::state_transition;
use koenig_planner::dynamics::AbsoluteOrbit;
use nalgebra::{Matrix6, SVector};
use std::f64::consts::PI;

fn wrap(x: f64) -> f64 {
    let mut y = x % (2.0 * PI);
    if y >= PI {
        y -= 2.0 * PI;
    }
    if y < -PI {
        y += 2.0 * PI;
    }
    y
}

/// eq.51 ROE map: x = [da, dlam, dex, dey, dix, diy] from chief & deputy mean elems.
fn roe(c: &AbsoluteOrbit, d: &AbsoluteOrbit) -> SVector<f64, 6> {
    let eta_c = (1.0 - c.e * c.e).sqrt();
    let dm = wrap(d.mean_anom - c.mean_anom);
    let dw = wrap(d.argp - c.argp);
    let dom = wrap(d.raan - c.raan);
    SVector::<f64, 6>::from_row_slice(&[
        (d.a - c.a) / c.a,
        dm + eta_c * (dw + dom * c.i.cos()),
        d.e * d.argp.cos() - c.e * c.argp.cos(),
        d.e * d.argp.sin() - c.e * c.argp.sin(),
        d.i - c.i,
        dom * c.i.sin(),
    ])
}

/// Pack/unpack AbsoluteOrbit <-> 6-vector [a,e,i,raan,argp,M].
fn to_vec(o: &AbsoluteOrbit) -> SVector<f64, 6> {
    SVector::<f64, 6>::from_row_slice(&[o.a, o.e, o.i, o.raan, o.argp, o.mean_anom])
}
fn from_vec(v: &SVector<f64, 6>) -> AbsoluteOrbit {
    AbsoluteOrbit::new(v[0], v[1], v[2], v[3], v[4], v[5])
}

/// Central-difference column with relative step `h_rel` for a vector function `f`.
fn cd_col(
    base: &SVector<f64, 6>,
    k: usize,
    h_rel: f64,
    f: &dyn Fn(&SVector<f64, 6>) -> SVector<f64, 6>,
) -> SVector<f64, 6> {
    let h = base[k].abs().max(1.0) * h_rel;
    let mut hp = *base;
    hp[k] += h;
    let mut hm = *base;
    hm[k] -= h;
    (f(&hp) - f(&hm)) / (2.0 * h)
}

/// Richardson-extrapolated Jacobian: combine steps h and h/2 to kill the O(h^2)
/// central-difference error, giving ~O(h^4) accuracy.
fn jacobian(
    base: &SVector<f64, 6>,
    f: &dyn Fn(&SVector<f64, 6>) -> SVector<f64, 6>,
) -> Matrix6<f64> {
    let mut j = Matrix6::zeros();
    for k in 0..6 {
        let c1 = cd_col(base, k, 1e-6, f);
        let c2 = cd_col(base, k, 5e-7, f);
        j.set_column(k, &((4.0 * c2 - c1) / 3.0));
    }
    j
}

/// d roe(chief, deputy) / d deputy at deputy = `dep`.
fn jac_roe(chief: &AbsoluteOrbit, dep: &AbsoluteOrbit) -> Matrix6<f64> {
    jacobian(&to_vec(dep), &|v| roe(chief, &from_vec(v)))
}

/// d propagate(oe, dt) / d oe at `oe`.
fn jac_propagate(oe: &AbsoluteOrbit, dt: f64) -> Matrix6<f64> {
    jacobian(&to_vec(oe), &|v| to_vec(&from_vec(v).propagate(dt)))
}

/// Independent FD reconstruction of Phi(t, t_f) about the chief (zero ROE).
fn fd_phi(chief_t: &AbsoluteOrbit, dt: f64) -> Matrix6<f64> {
    let chief_tf = chief_t.propagate(dt);
    let j_t = jac_roe(chief_t, chief_t); // d x(t)/d dep_t
    let j_f = jac_roe(&chief_tf, &chief_tf); // d x(t_f)/d dep_tf
    let p = jac_propagate(chief_t, dt); // d dep_tf / d dep_t
    j_f * p * j_t.try_inverse().expect("J_t invertible")
}

fn report(name: &str, chief: &AbsoluteOrbit, dt: f64) -> f64 {
    let analytic = state_transition(chief, &chief.propagate(dt), dt);
    let fd = fd_phi(chief, dt);
    // Noise floor: ignore entries far below the matrix scale, where FD roundoff
    // (~1e-9) dominates a near-zero analytic entry.
    let floor = 1e-7;
    let mut worst = 0.0f64;
    println!("=== {name} (dt={dt}) ===");
    for r in 0..6 {
        for c in 0..6 {
            let a = analytic[(r, c)];
            let f = fd[(r, c)];
            if a.abs().max(f.abs()) < floor {
                continue;
            }
            let rel = (a - f).abs() / a.abs().max(f.abs());
            if rel > 5e-5 {
                println!(
                    "  MISMATCH Phi[{}{}]: analytic={:+.6e} fd={:+.6e} rel={:.2e}",
                    r + 1,
                    c + 1,
                    a,
                    f,
                    rel
                );
            }
            worst = worst.max(rel);
        }
    }
    println!("  worst relative mismatch (above {floor:.0e} floor) = {worst:.3e}");
    worst
}

#[test]
fn stm_matches_independent_finite_difference() {
    // Worked-example chief at t=16050, dt = t_f - t.
    let chief0 = AbsoluteOrbit::new(
        25_000e3,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let w1 = report(
        "worked-example chief @ t=16050",
        &chief0.propagate(16_050.0),
        101_940.0,
    );
    let w2 = report("worked-example chief @ t=0", &chief0, 117_990.0);
    // Independent low-e fixture (the stm.rs oracle fixture).
    let fix = AbsoluteOrbit::new(
        25_000e3,
        0.3,
        50.0_f64.to_radians(),
        20.0_f64.to_radians(),
        40.0_f64.to_radians(),
        70.0_f64.to_radians(),
    );
    let w3 = report("low-e fixture", &fix, 39_000.0);

    // Hunter & D'Amico 2025 chief (e=0.70, i=51 deg, omega~200 deg): omega != 0
    // makes e_{y1} = e sin(omega) != 0, activating the Phi_34 / Phi_44 / Phi_64
    // delta-e couplings that are identically zero at the worked-example chief.
    let (hx, hy) = (-0.658_f64, -0.239_f64);
    let hunter = AbsoluteOrbit::new(
        25_000e3,
        (hx * hx + hy * hy).sqrt(),
        51.0_f64.to_radians(),
        30.0_f64.to_radians(),
        hy.atan2(hx),
        65.0_f64.to_radians() - hy.atan2(hx),
    );
    let w4 = report("Hunter chief @ t=0", &hunter, 39_000.0);
    let w5 = report(
        "Hunter chief @ t=10000",
        &hunter.propagate(10_000.0),
        29_000.0,
    );

    // FD precision floor on the smallest term (Phi_24, ~1e-5) is ~3e-5 through
    // the J_t^-1 amplification; a real coefficient/sign error is O(1e-2..1e0)
    // relative, so the 1e-4 threshold cleanly separates FD noise from a true mismatch.
    assert!(w1 < 1e-4, "worst @16050 = {w1:e}");
    assert!(w2 < 1e-4, "worst @0 = {w2:e}");
    assert!(w3 < 1e-4, "worst fixture = {w3:e}");
    assert!(w4 < 1e-4, "worst Hunter @0 = {w4:e}");
    assert!(w5 < 1e-4, "worst Hunter @10000 = {w5:e}");
}
