//! Adversarial-JSON tests over the central untrusted-string seam `run_json`.
//! Contract: any input yields Ok or a typed ApiError — never a panic, hang,
//! or OOM. `kind` is BadRequest or Solver (Internal is the HTTP panic-catch
//! kind and must not arise here).

use koenig_damico_planner_api::{run_json, ApiError, ApiErrorKind};
use proptest::prelude::*;

fn cost_spec_json() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(r#"{"type":"norm2"}"#.to_string()),
        Just(r#"{"type":"facemax"}"#.to_string()),
        Just(r#"{"type":"piecewise","period":null,"t_perigee0":null}"#.to_string()),
        Just(r#"{"type":"piecewise","period":39283.0,"t_perigee0":19641.0}"#.to_string()),
    ]
}

prop_compose! {
    // Request-shaped JSON spanning valid AND out-of-domain scalars. Spans/dt
    // are kept small so valid grids stay tiny (fast); extreme-dt DoS is in the
    // corpus (Task 7). Angles are degrees (the api wire convention).
    fn request_shaped_json()(
        a in -1.0e8..1.0e8f64,                 // includes a ≤ 0 (invalid)
        e in -0.5..1.5f64,                     // includes e ≥ 1 and e < 0 (invalid)
        i in -10.0..190.0f64,                  // includes sin i ≈ 0 (invalid)
        raan in -360.0..360.0f64,
        argp in -360.0..360.0f64,
        mean_anom in -360.0..360.0f64,
        t_i in 0.0..100.0f64,
        t_f in 0.0..100.0f64,
        dt in prop_oneof![-50.0..0.0f64, 0.5..50.0f64], // invalid or moderate
        w in proptest::array::uniform6(-1.0e3..1.0e3f64),
        cost in cost_spec_json(),
        init_len in 0usize..8,
    ) -> String {
        let w_str = w.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let init_str = (0..init_len)
            .map(|k| (t_i + k as f64).to_string())
            .collect::<Vec<_>>()
            .join(",");
        format!(
            r#"{{"chief":{{"a":{a},"e":{e},"i":{i},"raan":{raan},"argp":{argp},"mean_anom":{mean_anom}}},"t_i":{t_i},"t_f":{t_f},"dt":{dt},"w_metres":[{w_str}],"cost":{cost},"initial_times":[{init_str}]}}"#
        )
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn run_json_is_total(s in request_shaped_json()) {
        match run_json(&s) {
            Ok(_) => {}
            Err(ApiError { kind, .. }) => {
                prop_assert!(
                    matches!(kind, ApiErrorKind::BadRequest | ApiErrorKind::Solver),
                    "unexpected kind {:?} for input {}", kind, s
                );
            }
        }
    }
}
