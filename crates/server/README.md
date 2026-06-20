# koenig-damico-planner-server

Self-hostable HTTP service for the
[Koenig-D'Amico](https://github.com/sakobu/koenig-planner) fuel-optimal impulsive
maneuver planner (IEEE TAC 2020). A thin `axum` wrapper over
[`koenig-damico-planner-api`](../api) — the solver runs natively on your machine;
nothing is sent anywhere.

> Workspace-internal crate (`publish = false`): build it from this workspace, not
> from crates.io.

## Run

```bash
cargo run -p koenig-damico-planner-server
# listens on 0.0.0.0:8080 by default; override with KOENIG_PLANNER_ADDR.
# log level via RUST_LOG (default: info).
```

## Endpoints

| Method + path | Body                       | Success               |
| ------------- | -------------------------- | --------------------- |
| `GET /health` | —                          | `200 {"status":"ok"}` |
| `POST /solve` | `SolveRequest` (see below) | `200 SolveResponse`   |

The `SolveRequest` / `SolveResponse` JSON shapes are the
[`koenig-damico-planner-api`](../api) contract: chief angles in **degrees**,
`a` and `w_metres` in **metres**, times in **seconds**.

```bash
curl -s -H 'content-type: application/json' \
     -d @crates/server/golden.json \
     localhost:8080/solve | jq
```

## Errors

Every response is JSON. On failure the body is `{"kind": …, "message": …}`:

| `kind`        | Meaning                                              | Status                                                                 |
| ------------- | ---------------------------------------------------- | ---------------------------------------------------------------------- |
| `bad_request` | invalid input / malformed request body               | `400` (or the rejection's `415`/`422` for content-type / field errors) |
| `solver`      | well-formed request, numerically unsolvable / failed | `422`                                                                  |

The `kind` field is the source of truth: a `422` with `kind:"bad_request"` is a
request-field error from the extractor, whereas `kind:"solver"` is a planner
failure. Internal faults map to `500`. Requests are capped at 64 KiB; CORS is
permissive (the service is meant to be self-hosted, including by a browser demo).

## Docker

```bash
docker build -f crates/server/Dockerfile -t koenig-damico-planner-server .
docker run --rm -p 8080:8080 koenig-damico-planner-server
curl -s -d @crates/server/golden.json -H 'content-type: application/json' \
     localhost:8080/solve | jq
```

## License

Licensed under either of [Apache-2.0](../../LICENSE-APACHE) or [MIT](../../LICENSE-MIT)
at your option.
