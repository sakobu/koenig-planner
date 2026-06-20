"""Type stubs for the koenig_planner extension (PEP 561)."""

from typing import Sequence

__all__ = ["Orbit", "Maneuver", "Solution", "solve", "solve_json", "__version__"]

__version__: str

class Orbit:
    """Chief mean absolute orbit. Angles in degrees; `a` in metres."""

    a: float
    e: float
    i: float
    raan: float
    argp: float
    mean_anom: float
    def __init__(
        self,
        a: float,
        e: float,
        i: float,
        raan: float,
        argp: float,
        mean_anom: float,
    ) -> None: ...

class Maneuver:
    """One impulsive maneuver: time `t` [s] and RTN delta-v [m/s]."""

    # Returned by solve(); not constructed directly.
    t: float
    dv: tuple[float, float, float]

class Solution:
    """Planner output."""

    # Returned by solve(); not constructed directly.
    maneuvers: list[Maneuver]
    total_dv: float
    iterations: int
    residual: float
    lambda_: list[float]

def solve(
    chief: Orbit,
    t_i: float,
    t_f: float,
    dt: float,
    w_metres: Sequence[float],
    cost: str = ...,
    *,
    period: float | None = ...,
    t_perigee0: float | None = ...,
    n_coarse: int | None = ...,
    n_init: int | None = ...,
    eps_cost: float | None = ...,
    eps_remove: float | None = ...,
    initial_times: Sequence[float] | None = ...,
) -> Solution:
    """Plan a maneuver set. `cost` is one of "norm2", "facemax", "piecewise"."""
    ...

def solve_json(input: str) -> str:
    """Parse a JSON SolveRequest, run it, return the JSON SolveResponse."""
    ...
