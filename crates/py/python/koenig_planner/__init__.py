"""Python bindings for the Koenig-D'Amico fuel-optimal impulsive maneuver planner.

The public API is implemented in the compiled `_koenig_planner` extension and
re-exported here so `import koenig_planner` exposes it directly.
"""

from ._koenig_planner import (  # pyright: ignore[reportMissingModuleSource]  # native ext: stub shipped, binary built at install time
    Maneuver,
    Orbit,
    Solution,
    __version__,
    solve,
    solve_json,
)

__all__ = [
    "Orbit",
    "Maneuver",
    "Solution",
    "solve",
    "solve_json",
    "__version__",
]
