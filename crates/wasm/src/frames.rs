//! Coordinate-frame conversions and quasi-nonsingular ROE↔element algebra for
//! the presentation geometry. Reuses the core's `AbsoluteOrbit` (mean elements,
//! Kepler solver, J2 secular propagation); no dynamics are re-implemented. The
//! only math defined here is exact closed-form frame rotations and the exact
//! ROE definition/inverse, each pinned by the tests below.

use koenig_damico_planner_api::core::dynamics::AbsoluteOrbit;
use koenig_damico_planner_api::core::PlannerError;

// Private helpers are defined here for reuse in later geometry tasks; suppress
// dead_code until they are wired up.
#[allow(dead_code)]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[allow(dead_code)]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[allow(dead_code)]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[allow(dead_code)]
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Perifocal (PQW) → ECI (IJK) direction-cosine matrix from the orientation
/// angles `(i, Ω, ω)` `[rad]`, returned as the ECI column unit vectors
/// `[P̂, Q̂, Ŵ]` (perigee / semi-latus / orbit-normal). Standard 3-1-3 rotation
/// (e.g. Vallado, *Fundamentals of Astrodynamics*, perifocal-to-IJK DCM).
#[allow(dead_code)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
pub fn position_eci(orbit: &AbsoluteOrbit) -> Result<[f64; 3], PlannerError> {
    let nu = orbit.true_anomaly()?;
    Ok(orbit_point_eci(
        orbit.a, orbit.e, orbit.i, orbit.raan, orbit.argp, nu,
    ))
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
}
