import math
import pytest
import koenig_planner as kp

CHIEF = kp.Orbit(a=25_000e3, e=0.7, i=40.0, raan=358.0, argp=0.0, mean_anom=180.0)
W = [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0]
WINDOW = (0.0, 117_990.0, 30.0)


def test_golden_worked_example():
    sol = kp.solve(CHIEF, *WINDOW, W, "piecewise")
    assert 1 <= len(sol.maneuvers) <= 6
    assert 0.078 < sol.total_dv < 0.083
    assert sol.residual < 1e-3
    assert 1 <= sol.iterations <= 50
    assert len(sol.lambda_) == 6
    assert all(math.isfinite(x) for x in sol.lambda_)
    for m in sol.maneuvers:
        assert math.isfinite(m.t)
        assert len(m.dv) == 3 and all(math.isfinite(c) for c in m.dv)

    # Primer-vector history: three parallel, non-empty arrays; primer_rtn entries
    # are 3-tuples; magnitude reaches the |p| = 1 bound without exceeding
    # 1 + eps_cost (default 0.01), and is ~1 at each maneuver (complementary
    # slackness). This pins the PyO3 getters and the tuple shape the stub promises.
    n = len(sol.primer_times)
    assert n > 0
    assert len(sol.primer_magnitude) == n
    assert len(sol.primer_rtn) == n
    assert all(math.isfinite(g) for g in sol.primer_magnitude)
    assert all(
        len(p) == 3 and all(math.isfinite(c) for c in p) for p in sol.primer_rtn
    )
    max_g = max(sol.primer_magnitude)
    assert 1.0 - 1e-3 <= max_g <= 1.0 + 0.01 + 1e-6
    grid_times = sol.primer_times
    for m in sol.maneuvers:
        k = min(range(n), key=lambda i: abs(grid_times[i] - m.t))
        assert abs(grid_times[k] - m.t) < 1e-6
        assert 0.98 <= sol.primer_magnitude[k] <= 1.0 + 2e-2


def test_facemax_runs():
    sol = kp.solve(CHIEF, *WINDOW, W, "facemax")
    assert 1 <= len(sol.maneuvers) <= 6
    assert math.isfinite(sol.total_dv) and sol.total_dv > 0.0
    # Primer history is populated on the FaceMax path too (parallel arrays).
    assert len(sol.primer_times) > 0
    assert len(sol.primer_magnitude) == len(sol.primer_times)
    assert len(sol.primer_rtn) == len(sol.primer_times)
    assert all(len(p) == 3 for p in sol.primer_rtn)


def test_n_coarse_zero_is_value_error():
    with pytest.raises(ValueError):
        kp.solve(CHIEF, *WINDOW, W, "piecewise", n_coarse=0)


def test_unknown_cost_is_value_error():
    with pytest.raises(ValueError):
        kp.solve(CHIEF, *WINDOW, W, "bogus")


def test_orbit_field_access():
    assert CHIEF.a == 25_000e3 and CHIEF.e == 0.7 and CHIEF.mean_anom == 180.0
