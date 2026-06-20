# koenig-damico-planner-api

Internal serde / JSON facade for
[`koenig-damico-planner`](https://crates.io/crates/koenig-damico-planner) — the **one** place the
generic `solve` / `solve_from_initial_times` is monomorphized over the cost model. The Python
bindings (and the planned HTTP/WASM frontends) all call into it, so each frontend stays a thin,
plain-data wrapper.

> This is a **workspace-internal crate** (`publish = false`): depend on it by path within this
> workspace, not from crates.io. For the solver itself, use
> [`koenig-damico-planner`](https://crates.io/crates/koenig-damico-planner) directly.

## What it provides

- **`run(req: SolveRequest) -> Result<SolveResponse, ApiError>`** — the typed entry point.
- **`run_json(input: &str) -> Result<String, ApiError>`** — JSON-in / JSON-out convenience (parse a
  `SolveRequest`, run it, serialize the `SolveResponse`); used by the WASM/HTTP frontends.
- Plain-data DTOs with serde derives — `OrbitDto`, `CostSpec`, `SolveParamsDto`, `SolveRequest`,
  `ManeuverDto`, `SolveResponse`, `ApiError` — plus `pub use koenig_damico_planner as core;`.

`run()` owns the three unit/convention conversions so a frontend can't get them silently wrong:
request angles (`i`, `raan`, `argp`, `mean_anom`) are **degrees** → radians; the target `w_metres`
is divided by `chief.a`; the `Piecewise` cost period defaults to `TAU / chief.mean_motion()` when
omitted. Failures map to `ApiError { kind, message }` with `kind ∈ {"bad_request", "solver"}`.

## Usage

```rust,ignore
use koenig_damico_planner_api::{run, CostSpec, OrbitDto, SolveRequest};

let req = SolveRequest {
    chief: OrbitDto { a: 25_000e3, e: 0.7, i: 40.0, raan: 358.0, argp: 0.0, mean_anom: 180.0 },
    t_i: 0.0,
    t_f: 117_990.0,
    dt: 30.0,
    w_metres: [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0], // target pseudostate [m]
    cost: CostSpec::Piecewise { period: None, t_perigee0: None },
    params: None,
    initial_times: None,
};

let resp = run(req).expect("worked example should solve");
println!("{} maneuvers, total dv = {} m/s", resp.maneuvers.len(), resp.total_dv);
```

The same request as a JSON string can be driven through `run_json` instead, which is what the
Python `solve_json` and the planned WASM `solve` delegate to.

## License

Licensed under either of [Apache-2.0](../../LICENSE-APACHE) or [MIT](../../LICENSE-MIT) at your
option.
