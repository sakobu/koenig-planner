import json
import pytest
import koenig_planner as kp

REQ = {
    "chief": {"a": 25_000e3, "e": 0.7, "i": 40.0, "raan": 358.0, "argp": 0.0, "mean_anom": 180.0},
    "t_i": 0.0,
    "t_f": 117_990.0,
    "dt": 30.0,
    "w_meters": [50.0, 5000.0, 100.0, 100.0, 0.0, 400.0],
    "cost": {"type": "piecewise"},
}


def test_solve_json_roundtrip():
    out = json.loads(kp.solve_json(json.dumps(REQ)))
    assert 1 <= len(out["maneuvers"]) <= 6
    assert 0.078 < out["total_dv"] < 0.083
    assert out["residual"] < 1e-3
    assert len(out["lambda"]) == 6


def test_solve_json_malformed_is_value_error():
    with pytest.raises(ValueError):
        kp.solve_json("{ not json")
