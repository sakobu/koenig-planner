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
# request timeout (s) via KOENIG_PLANNER_TIMEOUT_SECS (default: 10).
# max simultaneous solves via KOENIG_PLANNER_MAX_CONCURRENCY (default: 64).
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
| `internal`    | unexpected internal fault — any handler/middleware panic, caught and logged server-side | `500`                       |

The `kind` field is the source of truth: a `422` with `kind:"bad_request"` is a
request-field error from the extractor, whereas `kind:"solver"` is a planner
failure. Internal faults map to `500`. CORS is permissive (the service is meant to
be self-hosted, including by a browser demo); see [Limits](#limits) for request
bounds.

## Limits

The service is hardened against unbounded-cost requests:

- **Grid size** is capped at 100 000 points; a larger `(t_f − t_i)/dt` is rejected
  with `400 {kind:"bad_request"}` *before* any solve allocation. Cost scales with
  the point count, not the request byte size.
- **Request timeout** — `KOENIG_PLANNER_TIMEOUT_SECS` (default 10) → `408`.
- **Concurrency limit** — `KOENIG_PLANNER_MAX_CONCURRENCY` (default 64) bounds
  simultaneous solves; excess requests queue.
- **Body size** is capped at 64 KiB → `413`.

The `408`/`413` transport rejections are plain status responses; the
`{kind,message}` JSON contract covers application-level errors (including the
grid-size `400`).

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
