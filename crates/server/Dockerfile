# syntax=docker/dockerfile:1

# ---- cargo-chef base: caches dependency compilation ----
FROM rust:1@sha256:44637ff22d0a6571a221bfaf137849711ad02ff4723dbb4736e297538f6a3e60 AS chef
RUN cargo install cargo-chef --locked
WORKDIR /app

# ---- plan: capture the dependency recipe ----
FROM chef AS planner
COPY . .
RUN cargo chef prepare --recipe-path recipe.json

# ---- build: cook deps (cached), then build the server ----
FROM chef AS builder
COPY --from=planner /app/recipe.json recipe.json
RUN cargo chef cook --release --locked -p koenig-damico-planner-server --recipe-path recipe.json
COPY . .
RUN cargo build --release --locked -p koenig-damico-planner-server

# ---- runtime: distroless nonroot (glibc; deps are pure-Rust) ----
FROM gcr.io/distroless/cc-debian12:nonroot@sha256:66aa873a4a14fb164aa01296058efd8253744606d72715e45acface073359faa AS runtime
COPY --from=builder /app/target/release/koenig-damico-planner-server /
# Run as the distroless nonroot user (UID 65532); port 8080 (>1024) needs no privilege.
USER nonroot:nonroot
EXPOSE 8080
ENV KOENIG_PLANNER_ADDR=0.0.0.0:8080
# No HEALTHCHECK: distroless has no shell. Probe GET /health from the orchestrator.
ENTRYPOINT ["/koenig-damico-planner-server"]
