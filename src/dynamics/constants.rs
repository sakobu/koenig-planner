//! Physical constants used by the J2 mean-ROE dynamics (\[KD20\] Appendix, eq. 50).

/// Earth gravitational parameter [m^3/s^2].
pub const MU: f64 = 3.986e14;

/// Earth equatorial radius `[m]`.
pub const R_E: f64 = 6.378e6;

/// J2 zonal harmonic coefficient (dimensionless).
pub const J2: f64 = 1.082e-3;

/// Shared numerator of the J2 secular prefactor, `3·J2·R_E²·√µ` `[m^(7/2)/s]`.
///
/// Both the secular angle rates (eq. 50) and the STM scaling `κ` (p. 13) start
/// from this chief-independent factor and then divide it by their own
/// `a^(7/2)`-and-`η` denominator. Extracting only the numerator — not the full
/// per-site prefactor — is deliberate: it leaves each call site's division
/// grouping untouched, so the f64 result stays bit-identical to the original
/// inline expression (the two sites' denominators are not the same quantity).
///
/// Ref: \[KD20\] eq. 50; p. 13 (STM scaling `κ`).
pub(crate) fn j2_secular_numerator() -> f64 {
    3.0 * J2 * R_E * R_E * MU.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    // The helper must reproduce, bit-for-bit, the inline numerator that
    // `orbit::secular_rates` and `stm::state_transition` previously spelled out,
    // so centralizing it cannot perturb the dynamics. Compare bit patterns (not
    // an approximate equality) because exact reproduction is the guarantee.
    #[test]
    fn j2_secular_numerator_is_bit_identical_to_inline() {
        assert_eq!(
            j2_secular_numerator().to_bits(),
            (3.0 * J2 * R_E * R_E * MU.sqrt()).to_bits()
        );
    }
}
