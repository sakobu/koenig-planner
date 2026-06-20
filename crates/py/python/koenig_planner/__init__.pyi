"""Type stubs for the koenig_planner package (PEP 561).

Re-exports the public API from the compiled `_koenig_planner` extension so
`import koenig_planner` is fully typed. The authoritative declarations live in
`_koenig_planner.pyi`.
"""

from ._koenig_planner import Maneuver as Maneuver
from ._koenig_planner import Orbit as Orbit
from ._koenig_planner import Solution as Solution
from ._koenig_planner import __version__ as __version__
from ._koenig_planner import solve as solve
from ._koenig_planner import solve_json as solve_json

__all__ = [
    "Orbit",
    "Maneuver",
    "Solution",
    "solve",
    "solve_json",
    "__version__",
]
