# %% [markdown]
# # Primer-vector history — Python showcase
#
# Plans the canonical worked example (Koenig & D'Amico 2020, Table III) entirely
# in Python via the `koenig_planner` extension, then plots the dual-gauge primer
# magnitude over the planning window (the paper's Fig. 7 contact curve) together
# with the |p| = 1 reference line and the optimal maneuver times. A second
# subplot shows the primer vector's RTN components. Open in Jupyter/VS Code
# (cells delimited by `# %%`) or run directly:
# `python crates/py/examples/primer_history.py`.

# %%
from pathlib import Path

import matplotlib  # pyright: ignore[reportMissingImports]  # optional viz dep, not a package requirement

matplotlib.use("Agg")  # headless-safe; remove for interactive use
import matplotlib.pyplot as plt  # pyright: ignore[reportMissingImports]

import koenig_planner as kp

# %%
chief = kp.Orbit(a=25_000e3, e=0.7, i=40.0, raan=358.0, argp=0.0, mean_anom=180.0)
w_metres = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]

sol = kp.solve(chief, 0.0, 117_990.0, 30.0, w_metres, "piecewise")

print(f"maneuvers    : {len(sol.maneuvers)}")
print(f"total_dv     : {sol.total_dv:.6f} m/s")
print(f"primer points: {len(sol.primer_times)}")
print(f"max |p|      : {max(sol.primer_magnitude):.6f} (dimensionless)")

# %%
times = sol.primer_times
mag = sol.primer_magnitude
maneuver_times = [m.t for m in sol.maneuvers]

r = [p[0] for p in sol.primer_rtn]
t = [p[1] for p in sol.primer_rtn]
n = [p[2] for p in sol.primer_rtn]

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)

# --- Primer magnitude (dimensionless) ---
ax0.plot(times, mag, color="tab:blue", label="|p(t)|")
ax0.axhline(1.0, color="tab:red", linestyle="--", linewidth=1.0, label="|p| = 1")
for k, mt in enumerate(maneuver_times):
    ax0.axvline(
        mt,
        color="tab:green",
        linestyle=":",
        linewidth=1.0,
        label="maneuver" if k == 0 else None,
    )
ax0.set_ylabel("primer magnitude [-]")  # dimensionless
ax0.legend(loc="lower right")
ax0.set_title("Primer-vector magnitude (Fig. 7 contact curve)")

# --- Primer vector RTN components ---
ax1.plot(times, r, label="R")
ax1.plot(times, t, label="T")
ax1.plot(times, n, label="N")
for mt in maneuver_times:
    ax1.axvline(mt, color="tab:green", linestyle=":", linewidth=1.0)
ax1.set_ylabel("primer p(t) RTN [-]")  # dimensionless
ax1.set_xlabel("time [s]")
ax1.legend(loc="lower right")
ax1.set_title("Primer vector RTN components")

fig.suptitle("Koenig-D'Amico primer-vector history (piecewise cost)")
fig.tight_layout()
out = Path(__file__).parent / "primer_history.png"
fig.savefig(out, dpi=120)
print(f"wrote {out}")
