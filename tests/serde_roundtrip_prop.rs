//! Property-based serde round-trip for the wire types (contract, not paper
//! math). JSON has no NaN/Inf, so generators use finite floats; serde_json is
//! lossless for finite f64, so equality is bit-exact.
#![cfg(feature = "serde")]

use koenig_damico_planner::dynamics::AbsoluteOrbit;
use koenig_damico_planner::{SolveParams, TimeGrid};
use proptest::prelude::*;

prop_compose! {
    fn any_orbit()(
        a in -1.0e12..1.0e12f64,
        e in -1.0e3..1.0e3f64,
        i in -1.0e3..1.0e3f64,
        raan in -1.0e3..1.0e3f64,
        argp in -1.0e3..1.0e3f64,
        mean_anom in -1.0e3..1.0e3f64,
    ) -> AbsoluteOrbit {
        AbsoluteOrbit::new(a, e, i, raan, argp, mean_anom)
    }
}

prop_compose! {
    fn any_grid()(
        t_i in -1.0e9..1.0e9f64,
        t_f in -1.0e9..1.0e9f64,
        dt in -1.0e6..1.0e6f64,
    ) -> TimeGrid {
        // Build the bare value type directly (uniform() would reject some of
        // these); round-trip must not depend on validity.
        TimeGrid { t_i, t_f, dt }
    }
}

prop_compose! {
    fn any_params()(
        n_coarse in 0usize..1000,
        n_init in 0usize..1000,
        eps_cost in -1.0e3..1.0e3f64,
        eps_remove in -1.0e3..1.0e3f64,
    ) -> SolveParams {
        SolveParams { n_coarse, n_init, eps_cost, eps_remove }
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 512, ..ProptestConfig::default() })]

    #[test]
    fn orbit_round_trips(o in any_orbit()) {
        let s = serde_json::to_string(&o).unwrap();
        let back: AbsoluteOrbit = serde_json::from_str(&s).unwrap();
        prop_assert_eq!(back, o);
    }

    #[test]
    fn grid_round_trips(g in any_grid()) {
        let s = serde_json::to_string(&g).unwrap();
        let back: TimeGrid = serde_json::from_str(&s).unwrap();
        prop_assert_eq!(back, g);
    }

    // SolveParams has no PartialEq, so assert serialize stability instead.
    #[test]
    fn params_serialize_is_stable(p in any_params()) {
        let s1 = serde_json::to_string(&p).unwrap();
        let back: SolveParams = serde_json::from_str(&s1).unwrap();
        let s2 = serde_json::to_string(&back).unwrap();
        prop_assert_eq!(s1, s2);
    }
}
