# koenig-planner (Python)

Python bindings for the [Koenig-D'Amico](https://github.com/sakobu/koenig-planner)
fuel-optimal impulsive maneuver planner (IEEE TAC 2020). The solver runs natively
(Rust) on your machine — nothing is sent anywhere.

## Install (from source)

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install maturin
maturin develop -m crates/py/Cargo.toml      # dev build into the venv
# or: maturin build --release -m crates/py/Cargo.toml   # build a wheel
```

## Usage

```python
import koenig_planner as kp

chief = kp.Orbit(a=25_000e3, e=0.7, i=40.0, raan=358.0, argp=0.0, mean_anom=180.0)
#   a [m]; i, raan, argp, mean_anom in DEGREES.

sol = kp.solve(
    chief,
    t_i=0.0, t_f=117_990.0, dt=30.0,        # planning window [s]
    w_meters=[50, 5000, 100, 100, 0, 400],  # target pseudostate [m]
    cost="piecewise",                        # "norm2" | "facemax" | "piecewise"
)

print(sol.total_dv, "m/s in", len(sol.maneuvers), "maneuvers")
for m in sol.maneuvers:
    print(m.t, m.dv)        # time [s], (R, T, N) [m/s]
print(sol.lambda_)          # optimal dual (6-vector)
print(sol.iterations, sol.residual)   # Algorithm-2 iterations; equality residual

# Primer-vector history (the paper's Fig. 7 contact curve) — three parallel
# arrays, one entry per grid point:
print(sol.primer_times[0], sol.primer_magnitude[0], sol.primer_rtn[0])
#   primer_rtn is the primer p(t) = Γᵀλ in RTN, NOT the executed thrust
#   direction — the two coincide only for the "norm2" cost.
```

`solve()` also takes optional keyword-only tuning arguments: `n_coarse` /
`n_init` (Algorithm-1 seeding), `eps_cost` / `eps_remove` (refinement
tolerances), the `"piecewise"`-only `period` / `t_perigee0`, and an
`initial_times` override (which bypasses Algorithm 1). See the type stub for the
full signature and defaults.

`solve_json(str) -> str` accepts/returns the JSON `SolveRequest`/`SolveResponse`
contract from `koenig-damico-planner-api`. Invalid input raises `ValueError`; solver
failures raise `RuntimeError`.

See `examples/showcase.py` for a plotting walkthrough, and
`examples/primer_history.py` for the primer-vector (Fig. 7 contact-curve) plot.

## Types & editor setup

The package ships PEP 561 type stubs (`py.typed` + `.pyi`), so editors and type
checkers (Pylance/pyright, mypy) get full autocomplete and checking for
`koenig_planner`. In VS Code, select the interpreter where you ran
`maturin develop` (Command Palette → "Python: Select Interpreter" → your `.venv`)
so imports resolve.
