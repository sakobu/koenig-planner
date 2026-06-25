"""Type stubs for the compiled koenig_planner extension module (PEP 561).

This mirrors the PyO3 surface defined in `crates/py/src/lib.rs`. Type checkers
read this stub; the real implementation is the compiled `_koenig_planner`
extension installed alongside it.

NOTE: unlike the Rust DTO boundaries, this stub is NOT compiler-enforced. When
a `#[pyclass]` field or `solve` signature changes in `crates/py/src/lib.rs`,
update the matching declaration here by hand — pyright checks usage, not parity.
"""

from typing import Sequence

__all__ = ["Orbit", "Maneuver", "Solution", "solve", "solve_json", "__version__"]

__version__: str

class Orbit:
    """Chief mean absolute orbit. Angles in degrees; `a` in metres."""

    def __init__(
        self,
        a: float,
        e: float,
        i: float,
        raan: float,
        argp: float,
        mean_anom: float,
    ) -> None: ...
    # Read-only (the underlying PyO3 fields are get-only).
    @property
    def a(self) -> float: ...
    @property
    def e(self) -> float: ...
    @property
    def i(self) -> float: ...
    @property
    def raan(self) -> float: ...
    @property
    def argp(self) -> float: ...
    @property
    def mean_anom(self) -> float: ...

class Maneuver:
    """One impulsive maneuver: time `t` [s] and RTN delta-v [m/s]."""

    # Returned by solve(); not constructed directly. All fields are read-only.
    @property
    def t(self) -> float: ...
    @property
    def dv(self) -> tuple[float, float, float]: ...

class Solution:
    """Planner output."""

    # Returned by solve(); not constructed directly. All fields are read-only.
    @property
    def maneuvers(self) -> list[Maneuver]: ...
    @property
    def total_dv(self) -> float:
        """Total fuel cost [m/s]: the minimized objective (the paper's "delta-v
        cost" c*) — sum of ||dv|| under the L2 cost, the polytope gauge sum(theta)
        under FaceMax."""
        ...
    @property
    def iterations(self) -> int: ...
    @property
    def residual(self) -> float: ...
    @property
    def lambda_(self) -> list[float]: ...
    # Primer-vector history (paper's Fig. 7 contact curve), parallel arrays, one
    # entry per grid point. Times in [s] from t_i; magnitude is dimensionless
    # (<= 1 + eps_cost, ~= 1 at maneuver times); primer_rtn is the primer vector
    # p(t) = Gamma^T(t)*lambda in RTN (not the executed thrust direction).
    @property
    def primer_times(self) -> list[float]: ...
    @property
    def primer_magnitude(self) -> list[float]: ...
    @property
    def primer_rtn(self) -> list[tuple[float, float, float]]: ...

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
    """Plan a maneuver set. `cost` is one of "norm2", "facemax", "piecewise".

    `n_coarse`/`n_init` are ignored when `initial_times` is supplied (that path
    bypasses Algorithm 1).
    """
    ...

def solve_json(input: str) -> str:
    """Parse a JSON SolveRequest, run it, return the JSON SolveResponse."""
    ...
