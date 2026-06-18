//! Independent, formula-free validation of `B(t)`.
//!
//! `control_input_matrix` implements the GVE control-input matrix from the
//! paper. This test reproduces the same Jacobian by a completely different
//! route that uses NONE of `B`'s formulas: convert the chief's mean elements to
//! Cartesian `r, v`; perturb the deputy's velocity by `±ε` along each RTN basis
//! vector; convert back to elements; form the eq.-51 ROE state; central-
//! difference. Agreement to ~1e-9 confirms `B` at the formula level (not merely
//! that the Rust matches the spec's transcription).

use koenig_planner::dynamics::b_matrix::control_input_matrix;
use koenig_planner::dynamics::AbsoluteOrbit;
use nalgebra::{Matrix3, SMatrix, SVector, Vector3};

const MU: f64 = 3.986e14;

fn wrap_pi(x: f64) -> f64 {
    use std::f64::consts::PI;
    (x + PI).rem_euclid(2.0 * PI) - PI
}

fn m_to_nu(m: f64, e: f64) -> f64 {
    let m = wrap_pi(m);
    let mut ecc = m + e * m.sin();
    for _ in 0..80 {
        let d = (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
        ecc -= d;
        if d.abs() < 1e-15 {
            break;
        }
    }
    ((1.0 - e * e).sqrt() * ecc.sin()).atan2(ecc.cos() - e)
}

fn nu_to_m(nu: f64, e: f64) -> f64 {
    let ecc = ((1.0 - e * e).sqrt() * nu.sin()).atan2(e + nu.cos());
    ecc - e * ecc.sin()
}

// 3-1-3 elementary rotations (perifocal -> ECI uses R3(-Om) R1(-i) R3(-w)).
fn r3(t: f64) -> Matrix3<f64> {
    let (c, s) = (t.cos(), t.sin());
    Matrix3::new(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0)
}
fn r1(t: f64) -> Matrix3<f64> {
    let (c, s) = (t.cos(), t.sin());
    Matrix3::new(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c)
}

fn coe_to_rv(a: f64, e: f64, i: f64, om: f64, w: f64, nu: f64) -> (Vector3<f64>, Vector3<f64>) {
    let p = a * (1.0 - e * e);
    let r = p / (1.0 + e * nu.cos());
    let r_pqw = Vector3::new(r * nu.cos(), r * nu.sin(), 0.0);
    let sp = (MU / p).sqrt();
    let v_pqw = Vector3::new(-sp * nu.sin(), sp * (e + nu.cos()), 0.0);
    let q = r3(-om) * r1(-i) * r3(-w);
    (q * r_pqw, q * v_pqw)
}

// Returns (a, e, i, Omega, omega, nu).
fn rv_to_coe(r: &Vector3<f64>, v: &Vector3<f64>) -> (f64, f64, f64, f64, f64, f64) {
    let rn = r.norm();
    let vn = v.norm();
    let h = r.cross(v);
    let hn = h.norm();
    let node = Vector3::new(0.0, 0.0, 1.0).cross(&h);
    let evec = ((vn * vn - MU / rn) * r - r.dot(v) * v) / MU;
    let e = evec.norm();
    let a = -MU / (2.0 * (vn * vn / 2.0 - MU / rn));
    let i = (h.z / hn).acos();
    let om = node.y.atan2(node.x);
    let w = (node.cross(&evec).dot(&h) / hn).atan2(node.dot(&evec));
    let nu = (evec.cross(r).dot(&h) / hn).atan2(evec.dot(r));
    (a, e, i, om, w, nu)
}

// eq.-51 ROE state of `dep` relative to `chief`, each = (a, e, i, Om, w, M).
fn roe(chief: &[f64; 6], dep: &[f64; 6]) -> SVector<f64, 6> {
    let (ac, ec, ic, omc, wc, mc) = (chief[0], chief[1], chief[2], chief[3], chief[4], chief[5]);
    let (ad, ed, id, omd, wd, md) = (dep[0], dep[1], dep[2], dep[3], dep[4], dep[5]);
    let etac = (1.0 - ec * ec).sqrt();
    SVector::<f64, 6>::from_row_slice(&[
        (ad - ac) / ac,
        wrap_pi(md - mc) + etac * (wrap_pi(wd - wc) + wrap_pi(omd - omc) * ic.cos()),
        ed * wd.cos() - ec * wc.cos(),
        ed * wd.sin() - ec * wc.sin(),
        id - ic,
        wrap_pi(omd - omc) * ic.sin(),
    ])
}

// B(t) reconstructed by central-differencing the ROE response to an RTN dv.
fn fd_b(orbit: &AbsoluteOrbit) -> SMatrix<f64, 6, 3> {
    let (a, e, i, om, w, m) = (
        orbit.a,
        orbit.e,
        orbit.i,
        orbit.raan,
        orbit.argp,
        orbit.mean_anom,
    );
    let nu_c = m_to_nu(m, e);
    let (rc, vc) = coe_to_rv(a, e, i, om, w, nu_c);
    let rhat = rc / rc.norm();
    let nh = rc.cross(&vc);
    let nhat = nh / nh.norm();
    let that = nhat.cross(&rhat); // T = N x R completes the right-handed triad
    let basis = [rhat, that, nhat]; // columns R, T, N
    let chief = [a, e, i, om, w, m];
    let eps = 1e-2; // m/s; central difference is ~1e-10 relative at this step

    let mut bfd = SMatrix::<f64, 6, 3>::zeros();
    for (j, b) in basis.iter().enumerate() {
        let (ap, ep, ip, omp, wp, nup) = rv_to_coe(&rc, &(vc + eps * b));
        let xp = roe(&chief, &[ap, ep, ip, omp, wp, nu_to_m(nup, ep)]);
        let (an, en, in_, omn, wn, nun) = rv_to_coe(&rc, &(vc - eps * b));
        let xm = roe(&chief, &[an, en, in_, omn, wn, nu_to_m(nun, en)]);
        let col = (xp - xm) / (2.0 * eps);
        for r in 0..6 {
            bfd[(r, j)] = col[r];
        }
    }
    bfd
}

fn frob_rel_err(got: &SMatrix<f64, 6, 3>, expected: &SMatrix<f64, 6, 3>) -> f64 {
    (got - expected).norm() / expected.norm()
}

#[test]
fn b_matrix_matches_independent_finite_difference() {
    // Well-conditioned fixture (e=0.3) and the worked-example chief at t=16050.
    let fixture = AbsoluteOrbit::new(
        25_000e3,
        0.3,
        50.0_f64.to_radians(),
        20.0_f64.to_radians(),
        40.0_f64.to_radians(),
        70.0_f64.to_radians(),
    );

    let chief = AbsoluteOrbit::new(
        25_000e3,
        0.7,
        40.0_f64.to_radians(),
        358.0_f64.to_radians(),
        0.0,
        180.0_f64.to_radians(),
    );
    let worked_at_16050 = chief.propagate(16_050.0);

    // Hunter & D'Amico 2025 chief (e=0.70, i=51 deg, omega~200 deg): exercises the
    // B(t) terms that vanish at the worked-example chief (omega=0 zeros sin(omega)).
    let (hx, hy) = (-0.658_f64, -0.239_f64);
    let hunter = AbsoluteOrbit::new(
        25_000e3,
        (hx * hx + hy * hy).sqrt(),
        51.0_f64.to_radians(),
        30.0_f64.to_radians(),
        hy.atan2(hx),
        65.0_f64.to_radians() - hy.atan2(hx),
    );

    for orbit in [fixture, worked_at_16050, hunter, hunter.propagate(20_000.0)] {
        let analytic = control_input_matrix(&orbit);
        let numeric = fd_b(&orbit);
        let rel = frob_rel_err(&numeric, &analytic);
        assert!(
            rel < 1e-6,
            "FD B disagrees with GVE B: Frobenius relative error {rel:.3e} (orbit a={}, e={})",
            orbit.a,
            orbit.e
        );
    }
}
