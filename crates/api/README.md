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

## Wire stability

The JSON request/response shape is part of this crate's public contract and is
versioned with the workspace crates: the wire schema only changes in a
semver-significant release. There is no `schema_version` field — the crate
version is the single source of truth for the contract.

Stable identifiers a client may hard-code:

- **Cost-model tags** (`cost.type`): `"norm2"`, `"facemax"`, `"piecewise"`.
- **Error kinds** (`ApiError.kind`): `"bad_request"`, `"solver"`, `"internal"`.
- **Field names** of `SolveRequest` / `SolveResponse` as documented on the DTOs.

These tags are regression-pinned by `crates/api/tests/serde_shapes.rs`, so a
silent rename can't slip through. Requests use `#[serde(deny_unknown_fields)]`:
an unknown request field is rejected rather than ignored. Responses may gain
fields in a future release, so clients should ignore unknown response fields.

## License

Licensed under either of [Apache-2.0](../../LICENSE-APACHE) or [MIT](../../LICENSE-MIT) at your
option.
