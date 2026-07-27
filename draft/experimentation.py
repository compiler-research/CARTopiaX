"""
Copyright 2026 compiler-research.org, Salvador de la Torre Gonzalez

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
    SPDX-License-Identifier: Apache-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This file contains the Bayesian optimization workflow used to calibrate
the model parameters, developed under Google Summer of Code (GSoC)
for the compiler-research.org organization.
"""

import json
import shutil
import subprocess
from pathlib import Path

# --- Adjust these paths/values to your environment ---
# Directory for the visualization configuration file bdm.toml
BDM_TOML_PATH = Path("./bdm.toml")
BASE_DIR = Path(__file__).resolve().parent.parent
# Directory for the parameters file params.json
PARAMS_PATH = BASE_DIR / "params.json"

BIODYNAMO_DIR = " /home/usuario/Desktop/biodynamo/build/bin/thisbdm.sh"
NUMBER_EXECUTIONS = 10

SIM_DIR = BASE_DIR / "output"        # folder where the ABM writes its results
DRAFT_DIR = BASE_DIR / "draft"       # base destination folder


def run_ABM(seed):

    config = {
        "seed": seed,
        "output_performance_statistics": False,
        "total_minutes_to_simulate": 4320,
        "output_csv_interval": 600,
        "bounded_space_length": 6500.0,
        "tumor_shape": "cylinder",
        "output_information_dependent_on_radius": True,
        "cylindrical_tumor_radius": 3000.0,
        "max_radius_analysis_csv_dependent_on_radius": 3000.0,
        "num_radius_intervals": 20,
        "cylindrical_tumor_height": 100.0,
        "initial_number_of_cylindrical_tumor_cells": 28000,
        "oncoprotein_mean": 1.0,
        "oncoprotein_standard_deviation": 0.0,
        "lateral_oxygen_production_min_z": -6600.0,
        "lateral_oxygen_production_max_z": 6600.0,
        "diffuse_on_z_axis": True,
        "initial_oxygen_level": 0.0,
        "oxygen_reference_level": 165.0,
        "default_oxygen_consumption_tumor_cell": 45.6,
        "diffusion_coefficient_oxygen": 180000.0,
        "decay_constant_oxygen": 0.01,
        "time_apoptosis": 6000.0,
        "treatment": {"0": 0},
        "average_time_transformation_random_rate": 72,
        "standard_deviation_transformation_random_rate": 15.0,
        "oxygen_saturation_for_proliferation": 0.0,
        "oxygen_limit_for_proliferation": 0.0,
        "oxygen_limit_for_necrosis_maximum": 0.0,
        "oxygen_limit_for_necrosis": 0.0,
        "reduction_consumption_dead_cells": 0.0,
        "basal_necrosis_probability_cancer_cells": 0.0000031,
        "nutrient_starvation_factor_cancer_cells": 1.000602,
        "dt_substances": 4320
    }

    # Save the config parameters for this run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)

    # Load the BioDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(
        ["bash", "-c", f"source {BIODYNAMO_DIR} && bdm run"], check=True
    )


def copy_output_to_draft(seed):
    """Copy the contents of SIM_DIR to draft/execution_seed_{seed}."""
    dest_dir = DRAFT_DIR / f"execution_seed_{seed}"

    if not SIM_DIR.exists():
        raise FileNotFoundError(f"Results folder not found: {SIM_DIR}")

    if dest_dir.exists():
        shutil.rmtree(dest_dir)  # avoid mixing results from previous runs

    shutil.copytree(SIM_DIR, dest_dir)
    print(f"[seed={seed}] Results copied to: {dest_dir}")


def run_multiple_seeds(N):
    """Run the ABM with seeds 0..N-1 and save each result separately."""
    DRAFT_DIR.mkdir(parents=True, exist_ok=True)

    # for seed in range(N):
    for seed in[4,5,6,7,8,9]:
        print(f"--- Running ABM with seed={seed} ---")
        run_ABM(seed)
        copy_output_to_draft(seed)

    print(f"Done: {N} runs completed (seeds 0 to {N - 1}).")

# Function to set the export value in the bdm.toml file
def set_export_value(path: Path, value: bool):
    lines = path.read_text().splitlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("export"):
            new_lines.append(f"export = {str(value).lower()}")
        else:
            new_lines.append(line)

    path.write_text("\n".join(new_lines))

if __name__ == "__main__":

    # Save the original content of the bdm.toml and params.json files to restore them later
    original_bdm = BDM_TOML_PATH.read_text()
    original_params = PARAMS_PATH.read_text()

    try:
        # Change the export value in the bdm.toml file to False to disable visualization
        set_export_value(BDM_TOML_PATH, False)

        run_multiple_seeds(N=NUMBER_EXECUTIONS)

    finally:
        # Restore the original content of the bdm.toml and params.json files always, even if an error occurs during the optimization process
        BDM_TOML_PATH.write_text(original_bdm)
        PARAMS_PATH.write_text(original_params)