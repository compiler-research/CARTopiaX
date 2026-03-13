"""Integration tests: run a short simulation end-to-end.

Requires BioDynaMo and CARTopiaX installed.  Set the environment variables
documented in cartopiaX/_bootstrap.py before running.

Run only these tests:
    pytest -m integration python/tests/test_integration.py
"""

import os

import pytest

pytestmark = pytest.mark.integration

if not os.environ.get("BDM_INSTALL_DIR"):
    pytest.skip("BDM_INSTALL_DIR not set", allow_module_level=True)

# ── Shared fixture ────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def sim_result(tmp_path_factory):
    """Run one short deterministic simulation; return its SimResult.

    Configuration choices:
      seed=42                   – reproducible RNG state
      initial_tumor_radius=50   – ~130 cells rather than ~3957; much faster
      treatment={}              – no CAR-T infusion; removes spawning code path
      total_minutes_to_simulate=720  – 12 hours; simulation runs 7201 steps
      output_csv_interval=1440  – write CSV every 144 minutes → ~5 data rows
      output_performance_statistics=False  – suppress per-step timing output
    """
    from cartopiaX import SimParam, Simulation

    run_dir = tmp_path_factory.mktemp("sim_run")
    out_dir = run_dir / "output"

    params = SimParam(
        seed=42,
        initial_tumor_radius=50.0,
        treatment={},
        total_minutes_to_simulate=720,
        output_csv_interval=1440,
        output_performance_statistics=False,
    )

    orig_dir = os.getcwd()
    os.chdir(run_dir)
    try:
        result = Simulation(params, output_dir=str(out_dir)).run()
    finally:
        os.chdir(orig_dir)

    return result


# ── Structural correctness tests ──────────────────────────────────────────────

def test_result_has_multiple_rows(sim_result):
    # With output_csv_interval=1440 and total_minutes=720 we expect rows at
    # steps 0 and 1440 at minimum — i.e. at least 2 rows.
    assert len(sim_result.csv_rows) >= 2, (
        f"Expected ≥2 CSV rows, got {len(sim_result.csv_rows)}"
    )


def test_all_expected_columns_present(sim_result):
    expected = {
        "total_days",
        "total_hours",
        "total_minutes",
        "tumor_radius",
        "num_cells",
        "num_tumor_cells",
        "tumor_cells_type1",
        "tumor_cells_type2",
        "tumor_cells_type3",
        "tumor_cells_type4",
        "tumor_cells_type5_dead",
        "num_alive_cart",
        "average_oncoprotein",
        "average_oxygen_cancer_cells",
    }
    actual = set(sim_result.csv_rows[0].keys())
    diff = expected.symmetric_difference(actual)
    assert not diff, f"Column mismatch: {diff}"


def test_survival_fraction_in_unit_interval(sim_result):
    sf = sim_result.survival_fraction
    assert 0.0 <= sf <= 1.0, f"survival_fraction={sf} is outside [0, 1]"


def test_no_cart_cells_without_treatment(sim_result):
    # treatment={} → SpawnCart never fires → num_alive_cart stays 0 throughout.
    for row in sim_result.csv_rows:
        cart = int(row["num_alive_cart"])
        assert cart == 0, (
            f"Found {cart} alive CAR-T cells at t={row['total_minutes']} min "
            "despite treatment={}"
        )


def test_first_row_is_time_zero(sim_result):
    assert float(sim_result.csv_rows[0]["total_minutes"]) == pytest.approx(0.0)


def test_time_is_monotonically_increasing(sim_result):
    minutes = [float(r["total_minutes"]) for r in sim_result.csv_rows]
    assert minutes == sorted(minutes), (
        f"total_minutes is not monotonically increasing: {minutes}"
    )


def test_final_counts_are_nonnegative(sim_result):
    assert sim_result.final_tumor_cells >= 0
    assert sim_result.final_alive_cart >= 0


def test_initial_tumor_cells_positive(sim_result):
    # A 50 µm radius sphere must contain at least one tumor cell.
    from cartopiaX.simulation import _alive_tumor
    first_alive = _alive_tumor(sim_result.csv_rows[0])
    assert first_alive > 0, "No alive tumor cells at t=0 — sphere may not have been seeded"


def test_tumor_radius_positive_at_t0(sim_result):
    r = float(sim_result.csv_rows[0]["tumor_radius"])
    assert r > 0.0, f"tumor_radius at t=0 is {r}, expected > 0"


def test_oxygen_level_in_physiological_range(sim_result):
    # Average oxygen across tumor cells should stay within 0–100 mmHg.
    for row in sim_result.csv_rows:
        o2 = float(row["average_oxygen_cancer_cells"])
        assert 0.0 <= o2 <= 100.0, (
            f"average_oxygen_cancer_cells={o2} mmHg at t={row['total_minutes']} "
            "is outside physiological range [0, 100]"
        )
