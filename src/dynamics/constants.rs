//! Physical constants used by the J2 mean-ROE dynamics (Appendix, eq. 50).

/// Earth gravitational parameter [m^3/s^2].
pub const MU: f64 = 3.986e14;

/// Earth equatorial radius [m].
pub const R_E: f64 = 6.378e6;

/// J2 zonal harmonic coefficient (dimensionless).
pub const J2: f64 = 1.082e-3;
