//! Coordinate-frame conversions and quasi-nonsingular ROE↔element algebra for
//! the presentation geometry. Reuses the core's `AbsoluteOrbit` (mean elements,
//! Kepler solver, J2 secular propagation); no dynamics are re-implemented. The
//! only math defined here is exact closed-form frame rotations and the exact
//! ROE definition/inverse, each pinned by the tests below.

use koenig_damico_planner_api::core::dynamics::kepler::wrap_to_pi;
use koenig_damico_planner_api::core::dynamics::AbsoluteOrbit;
use koenig_damico_planner_api::core::PlannerError;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Perifocal (PQW) → ECI (IJK) direction-cosine matrix from the orientation
/// angles `(i, Ω, ω)` `[rad]`, returned as the ECI column unit vectors
/// `[P̂, Q̂, Ŵ]` (perigee / semi-latus / orbit-normal). Standard 3-1-3 rotation
/// (e.g. Vallado, *Fundamentals of Astrodynamics*, perifocal-to-IJK DCM).
fn perifocal_to_eci(i: f64, raan: f64, argp: f64) -> [[f64; 3]; 3] {
    let (co, so) = (raan.cos(), raan.sin());
    let (cw, sw) = (argp.cos(), argp.sin());
    let (ci, si) = (i.cos(), i.sin());
    let p = [co * cw - so * sw * ci, so * cw + co * sw * ci, sw * si];
    let q = [-co * sw - so * cw * ci, -so * sw + co * cw * ci, cw * si];
    let w = [so * si, -co * si, ci];
    [p, q, w]
}

/// ECI position `[m]` on the ellipse at an explicit true anomaly `nu` `[rad]`,
/// for the instantaneous elements `(a, e, i, raan, argp)`. Orbit-*shape*
/// sampling: `r = a(1-e²)/(1 + e·cosν)`, `r_pqw = [r cosν, r sinν, 0]`, rotated
/// by `perifocal_to_eci`. No propagation.
pub fn orbit_point_eci(a: f64, e: f64, i: f64, raan: f64, argp: f64, nu: f64) -> [f64; 3] {
    let r = a * (1.0 - e * e) / (1.0 + e * nu.cos());
    let (xp, yp) = (r * nu.cos(), r * nu.sin());
    let [p, q, _] = perifocal_to_eci(i, raan, argp);
    [
        p[0] * xp + q[0] * yp,
        p[1] * xp + q[1] * yp,
        p[2] * xp + q[2] * yp,
    ]
}

/// ECI position `[m]` of the orbit at its current mean anomaly. Solves Kepler
/// via the core's `true_anomaly` (faithful by reuse), then `orbit_point_eci`.
///
/// # Errors
/// Propagates [`AbsoluteOrbit::true_anomaly`]'s error (non-elliptic `e`).
pub fn position_eci(orbit: &AbsoluteOrbit) -> Result<[f64; 3], PlannerError> {
    let nu = orbit.true_anomaly()?;
    Ok(orbit_point_eci(
        orbit.a, orbit.e, orbit.i, orbit.raan, orbit.argp, nu,
    ))
}

/// Orthonormal RTN basis at the orbit's current position, as ECI column
/// vectors `[R̂, T̂, N̂]`: `R̂` radial (position direction), `N̂` orbit-normal
/// (the `Ŵ` column of `perifocal_to_eci`, eccentricity-independent), and
/// `T̂ = N̂ × R̂` the **transverse** (cross-radial, along-track) axis —
/// perpendicular to `R̂` within the orbital plane and prograde. `T̂` is the
/// transverse, *not* the tangential/velocity direction: the two coincide only
/// at apsides and otherwise differ by the flight-path angle. This is the same
/// RTN frame the core's control-input matrix uses (the `B(t)` columns are R, T,
/// N; see `koenig_damico_planner::dynamics::b_matrix`), so RTN Δv and primer
/// vectors rotate into ECI consistently with the solution. Right-handed
/// (`R̂ × T̂ = N̂`).
///
/// # Errors
/// Propagates [`position_eci`]'s error.
pub fn rtn_basis_eci(orbit: &AbsoluteOrbit) -> Result<[[f64; 3]; 3], PlannerError> {
    let r_hat = normalize(position_eci(orbit)?);
    let [_, _, n_hat] = perifocal_to_eci(orbit.i, orbit.raan, orbit.argp);
    let t_hat = cross(n_hat, r_hat);
    Ok([r_hat, t_hat, n_hat])
}

/// Rotate an RTN vector `(R, T, N)` into ECI using the basis at `orbit`.
///
/// # Errors
/// Propagates [`rtn_basis_eci`]'s error.
pub fn rtn_to_eci(orbit: &AbsoluteOrbit, v_rtn: [f64; 3]) -> Result<[f64; 3], PlannerError> {
    let [r, t, n] = rtn_basis_eci(orbit)?;
    Ok([
        r[0] * v_rtn[0] + t[0] * v_rtn[1] + n[0] * v_rtn[2],
        r[1] * v_rtn[0] + t[1] * v_rtn[1] + n[1] * v_rtn[2],
        r[2] * v_rtn[0] + t[2] * v_rtn[1] + n[2] * v_rtn[2],
    ])
}

/// Project an ECI vector into the chief's RTN frame at `orbit` (transpose of
/// [`rtn_to_eci`]: components along `R̂, T̂, N̂`).
///
/// # Errors
/// Propagates [`rtn_basis_eci`]'s error.
pub fn eci_to_rtn(orbit: &AbsoluteOrbit, v_eci: [f64; 3]) -> Result<[f64; 3], PlannerError> {
    let [r, t, n] = rtn_basis_eci(orbit)?;
    Ok([dot(v_eci, r), dot(v_eci, t), dot(v_eci, n)])
}

/// Reconstruct the deputy's mean absolute orbit from the chief and a
/// dimensionless quasi-nonsingular ROE offset `[δa, δλ, δeₓ, δe_y, δiₓ, δi_y]`
/// (`[KD20]` eq. 51, with the modified relative mean longitude
/// `δλ = δM + η·(δω + δΩ·cos i)`, `η = √(1−e²)` — the same convention as the
/// core's Φ/B matrices and FD oracles). This is the exact algebraic inverse of
/// the ROE definition — not a linearization — so the resulting deputy
/// `propagate`s with the same core dynamics as the chief.
pub fn deputy_from_roe(chief: &AbsoluteOrbit, roe: [f64; 6]) -> AbsoluteOrbit {
    let [da, dl, dex, dey, dix, diy] = roe;
    let a_d = chief.a * (1.0 + da);
    let i_d = chief.i + dix;
    let raan_d = chief.raan + diy / chief.i.sin();
    let ex_d = chief.e * chief.argp.cos() + dex;
    let ey_d = chief.e * chief.argp.sin() + dey;
    let e_d = (ex_d * ex_d + ey_d * ey_d).sqrt();
    let argp_d = ey_d.atan2(ex_d);
    // δω is wrapped because atan2 re-branches argp_d while the chief's argp is
    // unbounded after propagation; δΩ = δi_y/sin i is small by construction.
    let eta = (1.0 - chief.e * chief.e).sqrt();
    let d_argp = wrap_to_pi(argp_d - chief.argp);
    let d_raan = raan_d - chief.raan;
    let mean_anom_d = chief.mean_anom + dl - eta * (d_argp + d_raan * chief.i.cos());
    AbsoluteOrbit::new(a_d, e_d, i_d, raan_d, argp_d, mean_anom_d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use wasm_bindgen_test::wasm_bindgen_test;

    const A: f64 = 25_000e3;
    const E: f64 = 0.7;

    #[wasm_bindgen_test]
    fn perifocal_to_eci_equatorial_is_identity() {
        let [p, q, w] = perifocal_to_eci(0.0, 0.0, 0.0);
        assert!((p[0] - 1.0).abs() < 1e-12 && p[1].abs() < 1e-12 && p[2].abs() < 1e-12);
        assert!(q[0].abs() < 1e-12 && (q[1] - 1.0).abs() < 1e-12 && q[2].abs() < 1e-12);
        assert!(w[0].abs() < 1e-12 && w[1].abs() < 1e-12 && (w[2] - 1.0).abs() < 1e-12);
    }

    #[wasm_bindgen_test]
    fn perigee_radius_at_nu_zero_equatorial() {
        // i=Ω=ω=0, ν=0 → perigee on +X at radius a(1-e).
        let r = orbit_point_eci(A, E, 0.0, 0.0, 0.0, 0.0);
        assert!((r[0] - A * (1.0 - E)).abs() < 1e-3, "got {r:?}");
        assert!(r[1].abs() < 1e-6 && r[2].abs() < 1e-6, "got {r:?}");
    }

    #[wasm_bindgen_test]
    fn argp_rotates_perigee_into_plane() {
        // ω=90°, i=Ω=0, ν=0 → perigee points along +Y at radius a(1-e).
        let r = orbit_point_eci(A, E, 0.0, 0.0, PI / 2.0, 0.0);
        assert!(r[0].abs() < 1e-6 && r[2].abs() < 1e-6, "got {r:?}");
        assert!((r[1] - A * (1.0 - E)).abs() < 1e-3, "got {r:?}");
    }

    #[wasm_bindgen_test]
    fn position_eci_agrees_with_orbit_point() {
        let o = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.0, 1.2);
        let nu = o.true_anomaly().unwrap();
        let from_orbit = position_eci(&o).unwrap();
        let from_point = orbit_point_eci(o.a, o.e, o.i, o.raan, o.argp, nu);
        for k in 0..3 {
            assert!((from_orbit[k] - from_point[k]).abs() < 1e-6, "k={k}");
        }
    }

    #[wasm_bindgen_test]
    fn radius_bounded_by_perigee_and_apogee() {
        let o = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.0, 1.2);
        let r = norm(position_eci(&o).unwrap());
        assert!(
            (A * (1.0 - E) - 1.0..=A * (1.0 + E) + 1.0).contains(&r),
            "r={r}"
        );
    }

    #[wasm_bindgen_test]
    fn rtn_basis_is_orthonormal_right_handed() {
        let o = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.3, 1.2);
        let [r, t, n] = rtn_basis_eci(&o).unwrap();
        for v in [r, t, n] {
            assert!((norm(v) - 1.0).abs() < 1e-12);
        }
        assert!(dot(r, t).abs() < 1e-12 && dot(t, n).abs() < 1e-12 && dot(n, r).abs() < 1e-12);
        // right-handed: R̂ × T̂ = N̂
        let rt = cross(r, t);
        for k in 0..3 {
            assert!((rt[k] - n[k]).abs() < 1e-12, "k={k}");
        }
    }

    #[wasm_bindgen_test]
    fn radial_axis_points_along_position() {
        let o = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.3, 1.2);
        let [r, _, _] = rtn_basis_eci(&o).unwrap();
        let pos_hat = normalize(position_eci(&o).unwrap());
        for k in 0..3 {
            assert!((r[k] - pos_hat[k]).abs() < 1e-12, "k={k}");
        }
    }

    #[wasm_bindgen_test]
    fn normal_axis_is_equatorial_z_for_zero_inclination() {
        let o = AbsoluteOrbit::new(A, E, 0.0, 0.0, 0.0, 1.2);
        let [_, _, n] = rtn_basis_eci(&o).unwrap();
        assert!(n[0].abs() < 1e-12 && n[1].abs() < 1e-12 && (n[2] - 1.0).abs() < 1e-12);
    }

    #[wasm_bindgen_test]
    fn rtn_to_eci_preserves_magnitude() {
        let o = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.3, 1.2);
        let v: [f64; 3] = [1.5, -2.5, 0.75];
        let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((norm(rtn_to_eci(&o, v).unwrap()) - mag).abs() < 1e-12);
    }

    #[wasm_bindgen_test]
    fn rtn_eci_round_trip() {
        let o = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.3, 1.2);
        let v: [f64; 3] = [1.5, -2.5, 0.75];
        let back = eci_to_rtn(&o, rtn_to_eci(&o, v).unwrap()).unwrap();
        for k in 0..3 {
            assert!((back[k] - v[k]).abs() < 1e-12, "k={k}");
        }
    }

    fn chief_fixture() -> AbsoluteOrbit {
        AbsoluteOrbit::new(
            A,
            E,
            40f64.to_radians(),
            358f64.to_radians(),
            0.0,
            180f64.to_radians(),
        )
    }

    // Forward quasi-nonsingular ROE (deputy − chief), [KD20] eq. 51 with the
    // modified relative mean longitude δλ = δM + η_c·(δω + δΩ·cos i_c) — the
    // same ground-truth map the core's FD oracles pin (tests/fd_stm.rs,
    // tests/fd_b_matrix.rs at the repo root).
    fn roe_of(chief: &AbsoluteOrbit, deputy: &AbsoluteOrbit) -> [f64; 6] {
        let eta = (1.0 - chief.e * chief.e).sqrt();
        let dm = wrap_to_pi(deputy.mean_anom - chief.mean_anom);
        let dw = wrap_to_pi(deputy.argp - chief.argp);
        let dom = wrap_to_pi(deputy.raan - chief.raan);
        [
            (deputy.a - chief.a) / chief.a,
            dm + eta * (dw + dom * chief.i.cos()),
            deputy.e * deputy.argp.cos() - chief.e * chief.argp.cos(),
            deputy.e * deputy.argp.sin() - chief.e * chief.argp.sin(),
            deputy.i - chief.i,
            dom * chief.i.sin(),
        ]
    }

    #[wasm_bindgen_test]
    fn zero_roe_reproduces_chief() {
        let c = chief_fixture();
        let d = deputy_from_roe(&c, [0.0; 6]);
        assert!((d.a - c.a).abs() < 1e-6);
        assert!((d.e - c.e).abs() < 1e-15);
        assert!((d.i - c.i).abs() < 1e-15);
        assert!((d.raan - c.raan).abs() < 1e-15);
        assert!((d.argp - c.argp).abs() < 1e-15);
        assert!((d.mean_anom - c.mean_anom).abs() < 1e-15);
    }

    #[wasm_bindgen_test]
    fn pure_delta_a_scales_semimajor_only() {
        let c = chief_fixture();
        let d = deputy_from_roe(&c, [1e-4, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((d.a - c.a * (1.0 + 1e-4)).abs() < 1e-6);
        assert!((d.i - c.i).abs() < 1e-15 && (d.raan - c.raan).abs() < 1e-15);
    }

    #[wasm_bindgen_test]
    fn deputy_from_roe_is_exact_inverse_of_roe_of() {
        // Build an arbitrary deputy, take its ROE, reconstruct — must round-trip.
        let c = chief_fixture();
        let deputy = AbsoluteOrbit::new(
            A * (1.0 + 2e-4),
            0.705,
            40.01f64.to_radians(),
            358.02f64.to_radians(),
            0.05,
            180.1f64.to_radians(),
        );
        let roe = roe_of(&c, &deputy);
        let back = deputy_from_roe(&c, roe);
        assert!((back.a - deputy.a).abs() < 1e-3, "a");
        assert!((back.e - deputy.e).abs() < 1e-12, "e");
        assert!((back.i - deputy.i).abs() < 1e-12, "i");
        assert!((back.raan - deputy.raan).abs() < 1e-12, "raan");
        assert!((back.argp - deputy.argp).abs() < 1e-12, "argp");
        assert!((back.mean_anom - deputy.mean_anom).abs() < 1e-12, "M");
    }

    #[wasm_bindgen_test]
    fn inversion_roundtrips_with_perigee_offset() {
        // δω ≠ 0 activates the η weighting: at e = 0.7 an unweighted inverse
        // misses δM by (1 − η)·δω ≈ 0.29·δω — ~10 orders above this tolerance.
        let c = AbsoluteOrbit::new(A, E, 40f64.to_radians(), 358f64.to_radians(), 0.6, 2.1);
        let d = AbsoluteOrbit::new(
            A * (1.0 + 2e-4),
            0.703,
            40.01f64.to_radians(),
            358.02f64.to_radians(),
            0.65,
            2.15,
        );
        let back = deputy_from_roe(&c, roe_of(&c, &d));
        assert!((back.a - d.a).abs() < 1e-3, "a");
        assert!((back.e - d.e).abs() < 1e-12, "e");
        assert!((back.i - d.i).abs() < 1e-12, "i");
        assert!((back.raan - d.raan).abs() < 1e-12, "raan");
        assert!((back.argp - d.argp).abs() < 1e-12, "argp");
        assert!((back.mean_anom - d.mean_anom).abs() < 1e-12, "M");
    }

    #[wasm_bindgen_test]
    fn inversion_roundtrips_across_the_argp_wrap() {
        // Chief perigee just below +π, deputy's just past −π (the same physical
        // direction): δω must be taken the short way around, not ±2π off.
        let c = AbsoluteOrbit::new(
            A,
            E,
            40f64.to_radians(),
            358f64.to_radians(),
            PI - 0.01,
            1.0,
        );
        let d = AbsoluteOrbit::new(
            A * (1.0 + 1e-4),
            0.701,
            40.005f64.to_radians(),
            358.01f64.to_radians(),
            -PI + 0.01,
            1.02,
        );
        let back = deputy_from_roe(&c, roe_of(&c, &d));
        assert!((back.argp - d.argp).abs() < 1e-12, "argp");
        assert!((back.mean_anom - d.mean_anom).abs() < 1e-12, "M");
        assert!((back.e - d.e).abs() < 1e-12, "e");
    }
}
