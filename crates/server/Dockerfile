# syntax=docker/dockerfile:1

# ---- cargo-chef base: caches dependency compilation ----
FROM rust:1 AS chef
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

# ---- runtime: distroless (glibc; deps are pure-Rust) ----
FROM gcr.io/distroless/cc-debian12 AS runtime
COPY --from=builder /app/target/release/koenig-damico-planner-server /
EXPOSE 8080
ENV KOENIG_PLANNER_ADDR=0.0.0.0:8080
ENTRYPOINT ["/koenig-damico-planner-server"]
