"""SimParam: Python representation of bdm::SimParam from hyperparams.h.

Ownership
---------
SimParam instances are owned by Python.  When the simulation runs, the params
are serialised to params.json; the C++ side constructs its own bdm::SimParam
via SimParam::LoadParams().  No C++ pointer is ever held by Python for this
type.

Treatment map
-------------
The C++ field is  std::map<int, int>  keyed by treatment day.  Here it is a
plain Python dict[int, int].  Serialisation uses string keys because
nlohmann::json represents JSON object keys as strings and C++ LoadParams calls
std::stoi() to recover the integer day index.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict


@dataclass
class SimParam:
    """All fields from bdm::SimParam; defaults match the C++ header defaults."""

    # ── General ──────────────────────────────────────────────────────────────
    seed: int = 1
    output_performance_statistics: bool = False
    total_minutes_to_simulate: int = 43200
    initial_tumor_radius: float = 150.0
    bounded_space_length: int = 1000

    # dict[day, n_cells].  Empty dict → no treatment.
    treatment: Dict[int, int] = field(default_factory=lambda: {0: 3957, 8: 3957})

    # ── Time steps ───────────────────────────────────────────────────────────
    dt_substances: float = 0.01
    dt_mechanics: float = 0.1
    dt_cycle: float = 6.0
    dt_step: float = 0.1
    # Steps between CSV writes.  C++ default: 12*60/dt_step = 7200 steps.
    output_csv_interval: int = 7200

    # ── Apoptosis volume change ───────────────────────────────────────────────
    volume_relaxation_rate_cytoplasm_apoptotic_cells: float = 0.0166667
    volume_relaxation_rate_nucleus_apoptotic_cells: float = 0.00583333
    volume_relaxation_rate_fluid_apoptotic_cells: float = 0.0
    time_apoptosis: float = 516.0
    reduction_consumption_dead_cells: float = 0.1

    # ── Chemicals ────────────────────────────────────────────────────────────
    resolution_grid_substances: int = 50
    diffusion_coefficient_oxygen: float = 100000.0
    decay_constant_oxygen: float = 0.1
    diffusion_coefficient_immunostimulatory_factor: float = 1000.0
    decay_constant_immunostimulatory_factor: float = 0.016
    oxygen_reference_level: float = 38.0
    initial_oxygen_level: float = 38.0
    oxygen_saturation: float = 30.0

    # ── Forces ───────────────────────────────────────────────────────────────
    cell_repulsion_between_tumor_tumor: float = 10.0
    cell_repulsion_between_cart_cart: float = 50.0
    cell_repulsion_between_cart_tumor: float = 50.0
    cell_repulsion_between_tumor_cart: float = 10.0
    max_relative_adhesion_distance: float = 1.25
    cell_adhesion_between_tumor_tumor: float = 0.4
    cell_adhesion_between_cart_cart: float = 0.0
    cell_adhesion_between_cart_tumor: float = 0.0
    cell_adhesion_between_tumor_cart: float = 0.0
    length_box_mechanics: int = 22
    dnew: float = 0.15
    dold: float = -0.05
    max_squared_distance_cart_moving_towards_tumor_cell: float = 317.746

    # ── Tumor cell ───────────────────────────────────────────────────────────
    rate_secretion_immunostimulatory_factor: float = 10.0
    saturation_density_immunostimulatory_factor: float = 1.0
    oncoprotein_mean: float = 1.0
    oncoprotein_standard_deviation: float = 0.25
    oxygen_saturation_for_proliferation: float = 38.0
    oxygen_limit_for_proliferation: float = 10.0
    oxygen_limit_for_necrosis: float = 5.0
    oxygen_limit_for_necrosis_maximum: float = 2.5
    time_lysis: float = 86400.0
    maximum_necrosis_rate: float = 0.00277778
    default_oxygen_consumption_tumor_cell: float = 10.0
    default_volume_new_tumor_cell: float = 2494.0
    default_volume_nucleus_tumor_cell: float = 540.0
    default_fraction_fluid_tumor_cell: float = 0.75
    average_time_transformation_random_rate: float = 38.6
    standard_deviation_transformation_random_rate: float = 3.7
    adhesion_time: float = 60.0
    oncoprotein_limit: float = 0.5
    oncoprotein_saturation: float = 2.0
    volume_relaxation_rate_alive_tumor_cell_cytoplasm: float = 0.00216667
    volume_relaxation_rate_alive_tumor_cell_nucleus: float = 0.00366667
    volume_relaxation_rate_alive_tumor_cell_fluid: float = 0.0216667
    volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell: float = 5.33333e-05
    volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell: float = 0.000216667
    volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell: float = 0.000833333
    volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell: float = 5.33333e-05
    volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell: float = 0.000216667
    volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell: float = 0.000833333
    threshold_cancer_cell_type1: float = 1.5
    threshold_cancer_cell_type2: float = 1.0
    threshold_cancer_cell_type3: float = 0.5
    threshold_cancer_cell_type4: float = 0.0

    # ── CAR-T cell ───────────────────────────────────────────────────────────
    # C++ default: dt_cycle * 10 * 24 * 60 = 6 * 10 * 24 * 60 = 86400 min.
    # The value in hyperparams.h (12342.86) is a legacy placeholder; LoadParams
    # recomputes it from dt_cycle when the key is absent from JSON.  We match
    # the computed default here so the JSON round-trip is stable.
    average_maximum_time_untill_apoptosis_cart: float = 86400.0
    default_oxygen_consumption_cart: float = 1.0
    default_volume_new_cart_cell: float = 2494.0
    kill_rate_cart: float = 0.06667
    adhesion_rate_cart: float = 0.013
    max_adhesion_distance_cart: float = 18.0
    min_adhesion_distance_cart: float = 14.0
    minimum_distance_from_tumor_to_spawn_cart: float = 50.0
    persistence_time_cart: float = 10.0
    migration_bias_cart: float = 0.5
    migration_speed_cart: float = 5.0
    elastic_constant_cart: float = 0.01

    # ── Serialisation ────────────────────────────────────────────────────────

    def to_json_dict(self) -> dict:
        """Return a dict ready for json.dump() that bdm::SimParam::LoadParams reads.

        Treatment keys are strings: nlohmann::json models JSON object keys as
        strings, and C++ LoadParams recovers integer days via std::stoi().
        """
        return {
            "seed": self.seed,
            "output_performance_statistics": self.output_performance_statistics,
            "total_minutes_to_simulate": self.total_minutes_to_simulate,
            "initial_tumor_radius": self.initial_tumor_radius,
            "bounded_space_length": self.bounded_space_length,
            "treatment": {str(day): dose for day, dose in self.treatment.items()},
            "dt_substances": self.dt_substances,
            "dt_mechanics": self.dt_mechanics,
            "dt_cycle": self.dt_cycle,
            "dt_step": self.dt_step,
            "output_csv_interval": self.output_csv_interval,
            "volume_relaxation_rate_cytoplasm_apoptotic_cells": self.volume_relaxation_rate_cytoplasm_apoptotic_cells,
            "volume_relaxation_rate_nucleus_apoptotic_cells": self.volume_relaxation_rate_nucleus_apoptotic_cells,
            "volume_relaxation_rate_fluid_apoptotic_cells": self.volume_relaxation_rate_fluid_apoptotic_cells,
            "time_apoptosis": self.time_apoptosis,
            "reduction_consumption_dead_cells": self.reduction_consumption_dead_cells,
            "resolution_grid_substances": self.resolution_grid_substances,
            "diffusion_coefficient_oxygen": self.diffusion_coefficient_oxygen,
            "decay_constant_oxygen": self.decay_constant_oxygen,
            "diffusion_coefficient_immunostimulatory_factor": self.diffusion_coefficient_immunostimulatory_factor,
            "decay_constant_immunostimulatory_factor": self.decay_constant_immunostimulatory_factor,
            "oxygen_reference_level": self.oxygen_reference_level,
            "initial_oxygen_level": self.initial_oxygen_level,
            "oxygen_saturation": self.oxygen_saturation,
            "cell_repulsion_between_tumor_tumor": self.cell_repulsion_between_tumor_tumor,
            "cell_repulsion_between_cart_cart": self.cell_repulsion_between_cart_cart,
            "cell_repulsion_between_cart_tumor": self.cell_repulsion_between_cart_tumor,
            "cell_repulsion_between_tumor_cart": self.cell_repulsion_between_tumor_cart,
            "max_relative_adhesion_distance": self.max_relative_adhesion_distance,
            "cell_adhesion_between_tumor_tumor": self.cell_adhesion_between_tumor_tumor,
            "cell_adhesion_between_cart_cart": self.cell_adhesion_between_cart_cart,
            "cell_adhesion_between_cart_tumor": self.cell_adhesion_between_cart_tumor,
            "cell_adhesion_between_tumor_cart": self.cell_adhesion_between_tumor_cart,
            "length_box_mechanics": self.length_box_mechanics,
            "dnew": self.dnew,
            "dold": self.dold,
            "max_squared_distance_cart_moving_towards_tumor_cell": self.max_squared_distance_cart_moving_towards_tumor_cell,
            "rate_secretion_immunostimulatory_factor": self.rate_secretion_immunostimulatory_factor,
            "saturation_density_immunostimulatory_factor": self.saturation_density_immunostimulatory_factor,
            "oncoprotein_mean": self.oncoprotein_mean,
            "oncoprotein_standard_deviation": self.oncoprotein_standard_deviation,
            "oxygen_saturation_for_proliferation": self.oxygen_saturation_for_proliferation,
            "oxygen_limit_for_proliferation": self.oxygen_limit_for_proliferation,
            "oxygen_limit_for_necrosis": self.oxygen_limit_for_necrosis,
            "oxygen_limit_for_necrosis_maximum": self.oxygen_limit_for_necrosis_maximum,
            "time_lysis": self.time_lysis,
            "maximum_necrosis_rate": self.maximum_necrosis_rate,
            "default_oxygen_consumption_tumor_cell": self.default_oxygen_consumption_tumor_cell,
            "default_volume_new_tumor_cell": self.default_volume_new_tumor_cell,
            "default_volume_nucleus_tumor_cell": self.default_volume_nucleus_tumor_cell,
            "default_fraction_fluid_tumor_cell": self.default_fraction_fluid_tumor_cell,
            "average_time_transformation_random_rate": self.average_time_transformation_random_rate,
            "standard_deviation_transformation_random_rate": self.standard_deviation_transformation_random_rate,
            "adhesion_time": self.adhesion_time,
            "oncoprotein_limit": self.oncoprotein_limit,
            "oncoprotein_saturation": self.oncoprotein_saturation,
            "volume_relaxation_rate_alive_tumor_cell_cytoplasm": self.volume_relaxation_rate_alive_tumor_cell_cytoplasm,
            "volume_relaxation_rate_alive_tumor_cell_nucleus": self.volume_relaxation_rate_alive_tumor_cell_nucleus,
            "volume_relaxation_rate_alive_tumor_cell_fluid": self.volume_relaxation_rate_alive_tumor_cell_fluid,
            "volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell": self.volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell,
            "volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell": self.volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell,
            "volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell": self.volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell,
            "volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell": self.volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell,
            "volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell": self.volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell,
            "volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell": self.volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell,
            "threshold_cancer_cell_type1": self.threshold_cancer_cell_type1,
            "threshold_cancer_cell_type2": self.threshold_cancer_cell_type2,
            "threshold_cancer_cell_type3": self.threshold_cancer_cell_type3,
            "threshold_cancer_cell_type4": self.threshold_cancer_cell_type4,
            "average_maximum_time_untill_apoptosis_cart": self.average_maximum_time_untill_apoptosis_cart,
            "default_oxygen_consumption_cart": self.default_oxygen_consumption_cart,
            "default_volume_new_cart_cell": self.default_volume_new_cart_cell,
            "kill_rate_cart": self.kill_rate_cart,
            "adhesion_rate_cart": self.adhesion_rate_cart,
            "max_adhesion_distance_cart": self.max_adhesion_distance_cart,
            "min_adhesion_distance_cart": self.min_adhesion_distance_cart,
            "minimum_distance_from_tumor_to_spawn_cart": self.minimum_distance_from_tumor_to_spawn_cart,
            "persistence_time_cart": self.persistence_time_cart,
            "migration_bias_cart": self.migration_bias_cart,
            "migration_speed_cart": self.migration_speed_cart,
            "elastic_constant_cart": self.elastic_constant_cart,
        }

    def write_json(self, path: str | os.PathLike = "params.json") -> None:
        """Write serialised params to `path` for bdm::SimParam::LoadParams."""
        with open(path, "w") as fh:
            json.dump(self.to_json_dict(), fh, indent=2)
