# %% [markdown]
# # Koenig-D'Amico maneuver planner — Python showcase
#
# Plans the canonical worked example (Koenig & D'Amico 2020, Table III) entirely
# in Python via the `koenig_planner` extension, then plots the per-maneuver RTN
# delta-v components and |dv| over the planning window. Open in Jupyter/VS Code
# (cells delimited by `# %%`) or run directly: `python crates/py/examples/showcase.py`.

# %%
import math
from pathlib import Path

import matplotlib  # pyright: ignore[reportMissingImports]  # optional viz dep, not a package requirement

matplotlib.use("Agg")  # headless-safe; remove for interactive use
import matplotlib.pyplot as plt  # pyright: ignore[reportMissingImports]

import koenig_planner as kp

# %%
chief = kp.Orbit(a=25_000e3, e=0.7, i=40.0, raan=358.0, argp=0.0, mean_anom=180.0)
w_metres = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]

sol = kp.solve(chief, 0.0, 117_990.0, 30.0, w_metres, "piecewise")

print(f"maneuvers : {len(sol.maneuvers)}")
print(f"total_dv  : {sol.total_dv:.6f} m/s")
print(f"iterations: {sol.iterations}")
print(f"residual  : {sol.residual:.3e}")
for k, m in enumerate(sol.maneuvers):
    mag = math.sqrt(sum(c * c for c in m.dv))
    print(f"  #{k}: t={m.t:8.1f} s  dv=({m.dv[0]:+.4e}, {m.dv[1]:+.4e}, {m.dv[2]:+.4e})  |dv|={mag:.4e}")

# %%
times = [m.t for m in sol.maneuvers]
r = [m.dv[0] for m in sol.maneuvers]
t = [m.dv[1] for m in sol.maneuvers]
n = [m.dv[2] for m in sol.maneuvers]
mag = [math.sqrt(sum(c * c for c in m.dv)) for m in sol.maneuvers]

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
width = max(300.0, 0.004 * (sol.maneuvers[-1].t - sol.maneuvers[0].t + 1.0))
ax0.bar([x - width for x in times], r, width=width, label="R")
ax0.bar(times, t, width=width, label="T")
ax0.bar([x + width for x in times], n, width=width, label="N")
ax0.set_ylabel("Δv component [m/s]")
ax0.legend()
ax0.set_title("Per-maneuver RTN Δv")

ax1.stem(times, mag)
ax1.set_ylabel("|Δv| [m/s]")
ax1.set_xlabel("time [s]")
ax1.set_title("Maneuver magnitudes over the planning window")

fig.tight_layout()
out = Path(__file__).parent / "koenig_plan.png"
fig.savefig(out, dpi=120)
print(f"wrote {out}")
