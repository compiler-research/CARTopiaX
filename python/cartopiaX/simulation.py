"""Simulation control wrapper, result dataclass, and CSV parser.

How Simulate() is called
------------------------
bdm::Simulate(int argc, const char** argv) reads params.json from the current
working directory, runs the full simulation, and writes output/final_data.csv.
We write params.json immediately before calling it, so the C++ side sees our
Python-configured parameters.

We pass argc=0, argv=nullptr.  BioDynaMo treats this identically to an
invocation with no command-line flags; all configuration comes from params.json
and the registered param group.

CSV format (from OutputSummary::operator() in utils_aux.cc)
------------------------------------------------------------
total_days, total_hours, total_minutes, tumor_radius, num_cells,
num_tumor_cells, tumor_cells_type1, tumor_cells_type2, tumor_cells_type3,
tumor_cells_type4, tumor_cells_type5_dead, num_alive_cart,
average_oncoprotein, average_oxygen_cancer_cells

The first row is written at simulation step 0 (t=0).  Subsequent rows are
written every output_csv_interval steps.

Survival fraction
-----------------
Defined as alive_tumor_last / alive_tumor_first where alive_tumor = sum of
type1 + type2 + type3 + type4 cells (type5 = dead, excluded).  Returns 0.0
when the initial alive count is 0 (degenerate: no tumor seeded).
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from typing import List

from .params import SimParam

_CSV_FILENAME = "final_data.csv"


@dataclass
class SimResult:
    """Structured output from a completed simulation."""

    # All rows from final_data.csv; each row is a dict keyed by column name.
    csv_rows: List[dict] = field(default_factory=list)

    # Alive tumor cells (types 1–4) in the last recorded timestep.
    final_tumor_cells: int = 0

    # Alive CAR-T cells in the last recorded timestep.
    final_alive_cart: int = 0

    # alive_tumor_last / alive_tumor_first; 0.0 if first row has 0 alive cells.
    survival_fraction: float = 0.0


def _alive_tumor(row: dict) -> int:
    return (
        int(row["tumor_cells_type1"])
        + int(row["tumor_cells_type2"])
        + int(row["tumor_cells_type3"])
        + int(row["tumor_cells_type4"])
    )


def _parse_csv(path: str) -> SimResult:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        return SimResult()

    first_alive = _alive_tumor(rows[0])
    last_alive = _alive_tumor(rows[-1])
    # Guard is reachable: first_alive == 0 when all t=0 cells are type0
    # (unclassified) or type5 (dead), which can happen with a degenerate seed.
    survival = last_alive / first_alive if first_alive > 0 else 0.0

    return SimResult(
        csv_rows=rows,
        final_tumor_cells=last_alive,
        final_alive_cart=int(rows[-1]["num_alive_cart"]),
        survival_fraction=survival,
    )


class Simulation:
    """Configure and run one CARTopiaX simulation.

    Example::

        from cartopiaX import SimParam, Simulation

        p = SimParam(seed=42, total_minutes_to_simulate=720, treatment={})
        result = Simulation(p).run()

    The optional `behavior_names` list accepts names registered with
    @register_behavior.  In this PoC the names are stored and accessible
    (e.g., for logging or future wiring); the integration path for passing
    Python behaviors into ongoing C++ agent loops requires a separate
    BioDynaMo operation that calls back into Python and is out of scope here.
    """

    def __init__(
        self,
        params: SimParam,
        behavior_names: list[str] | None = None,
        params_path: str | os.PathLike = "params.json",
        output_dir: str | os.PathLike = "output",
    ) -> None:
        self._params = params
        self._behaviors = list(behavior_names or [])
        self._params_path = str(params_path)
        self._output_dir = str(output_dir)

    @property
    def behaviors(self) -> list[str]:
        """Names of registered behaviors attached to this simulation."""
        return list(self._behaviors)

    def run(self) -> SimResult:
        """Write params.json, invoke bdm::Simulate(), parse and return results."""
        os.makedirs(self._output_dir, exist_ok=True)
        self._params.write_json(self._params_path)
        self._call_simulate()

        csv_path = os.path.join(self._output_dir, _CSV_FILENAME)
        if not os.path.exists(csv_path):
            raise FileNotFoundError(
                f"Expected simulation output at '{csv_path}' but the file was "
                "not created.  Check that bdm::Simulate() ran without error and "
                "that the output directory is writable."
            )
        return _parse_csv(csv_path)

    def _call_simulate(self) -> None:
        from . import _bootstrap  # noqa: PLC0415
        import cppyy  # noqa: PLC0415

        _bootstrap.require()

        # argc=0, argv=nullptr is intentionally safe here: BioDynaMo's Simulate()
        # only iterates over argv[1..argc-1] for command-line flag parsing, so
        # with argc=0 it never dereferences argv.  All simulation configuration
        # is read from params.json (written above in run()), not from argv.
        ret = cppyy.gbl.bdm.Simulate(0, cppyy.nullptr)
        if ret != 0:
            raise RuntimeError(
                f"bdm::Simulate() returned non-zero exit code {ret}"
            )
