"""Unit tests for SimParam.  No BioDynaMo required."""

import json

import pytest

from cartopiaX.params import SimParam


def test_defaults_match_cpp():
    p = SimParam()
    assert p.seed == 1
    assert p.output_performance_statistics is False
    assert p.total_minutes_to_simulate == 43200
    assert p.initial_tumor_radius == 150.0
    assert p.bounded_space_length == 1000
    assert p.treatment == {0: 3957, 8: 3957}
    assert p.dt_substances == 0.01
    assert p.dt_mechanics == 0.1
    assert p.dt_cycle == 6.0
    assert p.dt_step == 0.1
    assert p.output_csv_interval == 7200
    assert p.oxygen_limit_for_necrosis == 5.0
    assert p.kill_rate_cart == pytest.approx(0.06667)
    assert p.migration_bias_cart == 0.5


def test_field_assignment():
    p = SimParam()
    p.seed = 42
    p.total_minutes_to_simulate = 720
    p.initial_tumor_radius = 50.0
    p.treatment = {}
    assert p.seed == 42
    assert p.total_minutes_to_simulate == 720
    assert p.initial_tumor_radius == 50.0
    assert p.treatment == {}


def test_treatment_keys_are_strings_in_json():
    # nlohmann::json parses object keys as strings; LoadParams calls std::stoi.
    p = SimParam(treatment={0: 100, 14: 200, 28: 300})
    d = p.to_json_dict()
    assert d["treatment"] == {"0": 100, "14": 200, "28": 300}
    for key in d["treatment"]:
        assert isinstance(key, str), f"treatment key {key!r} must be str, got {type(key)}"


def test_empty_treatment_serialises():
    p = SimParam(treatment={})
    d = p.to_json_dict()
    assert d["treatment"] == {}


def test_json_round_trip(tmp_path):
    p = SimParam(seed=7, initial_tumor_radius=80.0, treatment={3: 500})
    path = tmp_path / "params.json"
    p.write_json(path)
    with open(path) as fh:
        loaded = json.load(fh)
    assert loaded["seed"] == 7
    assert loaded["initial_tumor_radius"] == 80.0
    assert loaded["treatment"] == {"3": 500}


def test_all_required_cpp_keys_present():
    # Verify every key that C++ LoadParams explicitly reads is in to_json_dict().
    required = {
        "seed",
        "output_performance_statistics",
        "total_minutes_to_simulate",
        "initial_tumor_radius",
        "bounded_space_length",
        "treatment",
        "dt_substances",
        "dt_mechanics",
        "dt_cycle",
        "dt_step",
        "output_csv_interval",
        "diffusion_coefficient_oxygen",
        "decay_constant_oxygen",
        "oxygen_reference_level",
        "initial_oxygen_level",
        "oncoprotein_mean",
        "oncoprotein_standard_deviation",
        "oxygen_limit_for_necrosis",
        "kill_rate_cart",
        "adhesion_rate_cart",
        "migration_bias_cart",
        "migration_speed_cart",
        "average_maximum_time_untill_apoptosis_cart",
    }
    d = SimParam().to_json_dict()
    missing = required - d.keys()
    assert not missing, f"to_json_dict() is missing keys: {sorted(missing)}"


def test_write_json_creates_file(tmp_path):
    p = SimParam()
    dest = tmp_path / "out.json"
    p.write_json(dest)
    assert dest.exists()
    content = json.loads(dest.read_text())
    assert content["seed"] == 1


def test_independent_instances_do_not_share_treatment():
    # dataclass default_factory must create a new dict per instance.
    a = SimParam()
    b = SimParam()
    a.treatment[99] = 1
    assert 99 not in b.treatment
