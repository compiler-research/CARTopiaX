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
import logging
import subprocess
from pathlib import Path

import numpy as np
import optuna
import pandas as pd

# Change this: File Parameters for the desired experiment
EXPERIMENT_ID = 7
SEED = 42
# You can set the number of trials to 0 to skip the optimization and just load the best result from the database
NUMBER_OF_TRIALS = 5000
# Number of Monte Carlo simulations to run for each trial. Use one for an aproximation of the error with a single montecarlo run
NUMBER_MONTE_CARLO = 1
# BioDynaMo directory to execute the comand source thisbdm.sh, you can change it to your own path
BIODYNAMO_DIR = " /home/usuario/Desktop/biodynamo/build/bin/thisbdm.sh"

# Other Hyperparameters for the optimization process
# LOGGING
logging.basicConfig(level=logging.INFO)

# Directory for the visualization configuration file bdm.toml
BDM_TOML_PATH = Path("./bdm.toml")
# Directory for the parameters file params.json
PARAMS_PATH = Path(__file__).resolve().parent.parent / "params.json"

# Directory for this experiment
EXPERIMENT_DIR = Path("abm_calibration") / f"experiment_{EXPERIMENT_ID}"
EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)


# Function to run the ABM with the given parameters
def run_ABM(params, seed):
    # Change this: parameter to be optimized in the ABM simulation
    basal_necrosis_probability_cancer_cells = params[
        "basal_necrosis_probability_cancer_cells"
    ]

    # Change this configuration for the ABM run
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
        "initial_number_of_cylindrical_tumor_cells": 2800,
        "oncoprotein_mean": 1.0,
        "oncoprotein_standard_deviation": 0.0,
        "lateral_oxygen_production_min_z": -300.0,
        "lateral_oxygen_production_max_z": 300.0,
        "diffuse_on_z_axis": False,
        "initial_oxygen_level": 0.0,
        "oxygen_reference_level": 165.0,
        "default_oxygen_consumption_tumor_cell": 45.6,
        "diffusion_coefficient_oxygen": 180000.0,
        "decay_constant_oxygen": 0.01,
        "time_apoptosis": 6000.0,
        "treatment": {
            "0": 0
        },
        "average_time_transformation_random_rate": 42.0,
        "standard_deviation_transformation_random_rate": 8.0,
        "oxygen_saturation_for_proliferation": 0.0,
        "oxygen_limit_for_proliferation": 0.0,
        "oxygen_limit_for_necrosis_maximum": 0.0,
        "oxygen_limit_for_necrosis": 0.0,
        "reduction_consumption_dead_cells": 0.0,
        "basal_necrosis_probability_cancer_cells": basal_necrosis_probability_cancer_cells
        }
    # Save the config parameters for the run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)

    # Load the ByoDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(["bash", "-c", f"source {BIODYNAMO_DIR} && bdm run"], check=True)


# Change This: Function to compute the error between the ABM simulation results and the experimental data
def compute_error():
#     sim_dir = Path(__file__).resolve().parent.parent / "output" / "data_dependent_on_radius_tumor.csv"

#     if not sim_dir.exists():
#         logging.error("Missing simulation CSV: %s", sim_dir)
#         return float("inf")

#     df_s = pd.read_csv(
#         sim_dir,
#         usecols=[
#             "total_minutes",
#             "average_oxygen_all_cells_radius_2850_to_3000",
#             "average_oxygen_all_cells_radius_0_to_150",
#         ],
#     )

#     # See the value at the minute 30 for the average oxygen level in the simulation data
#     row = df_s[df_s["total_minutes"] == 30].iloc[0]
#     value_border = row["average_oxygen_all_cells_radius_2850_to_3000"]
#     value_border_in_mol_m3 = value_border / 585  # Convert from mmHg to mol/m3
#     target_value_border = target_outer
#     value_center = row["average_oxygen_all_cells_radius_0_to_150"]
#     value_center_in_mol_m3 = value_center / 585  # Convert from mmHg to mol/m3
#     target_value_center = target_inner

#    # Debug print de los 4 valores 
#     print(f"Value border: {value_border_in_mol_m3}")
#     print(f"Target value border: {target_value_border}")
#     print(f"Value center: {value_center_in_mol_m3}")
#     print(f"Target value center: {target_value_center}")

#    # MSE
#     mse = ((value_border_in_mol_m3 - target_value_border) ** 2 + (value_center_in_mol_m3 - target_value_center) ** 2) / 2 
    
#     return float(mse)
    sim_dir = Path(__file__).resolve().parent.parent / "output" / "final_data.csv"


    if not sim_dir.exists():
        logging.error("Missing simulation CSV: %s", sim_dir)
        return float("inf")

    df_s = pd.read_csv(
        sim_dir, usecols=["total_minutes", "tumor_cells_type5_dead"]
    )

   # take the value at the minute 4320 for the number of dead tumor cells in the simulation data
    row = df_s[df_s["total_minutes"] == 4320].iloc[0]
    value_dead_cells = row["tumor_cells_type5_dead"]
    target_value_dead_cells = 310

    # Compute the mean squared error (MSE) between the number of dead tumor cells in the target and simulation data
    mse = (value_dead_cells - target_value_dead_cells) ** 2

    return float(mse)


    # dafault 
    # target = Path(__file__).resolve().parent / "target_data" / "final_data.csv"
    # sim_dir = Path(__file__).resolve().parent.parent / "output" / "final_data.csv"

    # if not target.exists():
    #     logging.error("Missing target CSV: %s", target)
    #     return float("inf")

    # if not sim_dir.exists():
    #     logging.error("Missing simulation CSV: %s", sim_dir)
    #     return float("inf")

    # df_t = pd.read_csv(target, usecols=["total_minutes", "average_oxygen_cancer_cells"])
    # df_s = pd.read_csv(
    #     sim_dir, usecols=["total_minutes", "average_oxygen_cancer_cells"]
    # )

    # # Merge the target and simulation dataframes on the "total_minutes" column to align the data for error computation
    # merged = pd.merge(
    #     df_t, df_s, on="total_minutes", how="inner", suffixes=("_t", "_s")
    # )

    # if merged.empty:
    #     logging.error("No common minutes between target and simulation")
    #     return float("inf")

    # # Compute the mean squared error (MSE) between the average oxygen levels in the target and simulation data
    # mse = (
    #     (
    #         merged["average_oxygen_cancer_cells_s"]
    #         - merged["average_oxygen_cancer_cells_t"]
    #     )
    #     ** 2
    # ).mean()
    # return float(mse)

# Function to run the ABM with the given parameters
def run_ABM2(params, seed, num_cells):
    # Change this: parameter to be optimized in the ABM simulation
    oxygen_reference_level = params[
        "oxygen_reference_level"
    ]
    default_oxygen_consumption_tumor_cell = params[
        "default_oxygen_consumption_tumor_cell"
    ]

    # Change this configuration for the ABM run
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
        "initial_number_of_cylindrical_tumor_cells": num_cells,
        "oncoprotein_standard_deviation": 0.0,
        "lateral_oxygen_production_min_z": -6600.0,
        "lateral_oxygen_production_max_z": 6600.0,
        "diffuse_on_z_axis": True,
        "initial_oxygen_level": 0.0,
        "oxygen_reference_level": oxygen_reference_level,
        "diffusion_coefficient_oxygen": 180000.0,
        "decay_constant_oxygen": 0.01,
        "treatment": {
            "0": 0
        },
        "default_oxygen_consumption_tumor_cell": default_oxygen_consumption_tumor_cell,
        "oxygen_limit_for_necrosis_maximum": 0.0,
        "oxygen_limit_for_necrosis": 0.0
        }

    # Save the config parameters for the run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)

    # Load the ByoDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(["bash", "-c", f"source {BIODYNAMO_DIR} && bdm run"], check=True)

def compute_error2(target_inner):
    sim_dir = Path(__file__).resolve().parent.parent / "output" / "data_dependent_on_radius_tumor.csv"

    if not sim_dir.exists():
        logging.error("Missing simulation CSV: %s", sim_dir)
        return float("inf")

    df_s = pd.read_csv(
        sim_dir,
        usecols=[
            "total_minutes",
            "average_oxygen_all_cells_radius_2850_to_3000",
            "average_oxygen_all_cells_radius_0_to_150",
        ],
    )

    # See the value at the minute 30 for the average oxygen level in the simulation data
    row = df_s[df_s["total_minutes"] == 30].iloc[0]
    value_center = row["average_oxygen_all_cells_radius_0_to_150"]
    value_center_in_mol_m3 = value_center / 585  # Convert from mmHg to mol/m3
    target_value_center = target_inner

   # Debug print de los 4 valores 
    print(f"Value center: {value_center_in_mol_m3}")
    print(f"Target value center: {target_value_center}")

    mse=0.
   # MSE
    if value_center_in_mol_m3 < target_value_center:
        mse = (value_center_in_mol_m3 - target_value_center) ** 2
    
    
    return float(mse)

# Objective function for the Optuna optimization process
def objective(trial):
    # Change this: Define the parameters to be optimized and their ranges
    params = {
        "basal_necrosis_probability_cancer_cells": trial.suggest_float("basal_necrosis_probability_cancer_cells", 0.0000112, 0.0000195, step=0.0000001),
    }

    logging.info(f"Trial {trial.number} | params={params}")

    # Compute the error as the average of the errors from multiple Monte Carlo simulations varying the seed
    total_error = 0
    for seed in np.random.randint(0, 10000, NUMBER_MONTE_CARLO):
        # print(f"Running ABM with seed {seed} no cup") 
        # run_ABM2(params, int(seed), num_cells=28000)
        # error = compute_error2(target_inner=0.09)
        # total_error += error
        # print(f"Error for seed {seed}: {error}")
        # print(f"Running ABM with seed {seed} and num_cells=28000 with cup")
        # run_ABM(params, int(seed), num_cells=28000)
        # error = compute_error(target_inner=0.009, target_outer=0.16)
        # total_error += error
        # print(f"Error for seed {seed}: {error}")
        print(f"Running ABM with seed {seed} and num_cells=2800 with cup")
        run_ABM(params, int(seed))
        error = compute_error()
        total_error += error
        print(f"Error for seed {seed}: {error}")

    error = total_error / NUMBER_MONTE_CARLO

    logging.info(f"Trial {trial.number} | error={error}")

    return error


# Fix the random seed for reproducibility
np.random.seed(SEED)

# Create an Optuna study to optimize the parameters of the ABM
study = optuna.create_study(
    study_name="abm_calibration",
    storage=f"sqlite:///{EXPERIMENT_DIR / 'abm_optuna.db'}",
    load_if_exists=True,
    direction="minimize",
    sampler=optuna.samplers.TPESampler(seed=SEED),
)


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

        study.optimize(objective, n_trials=NUMBER_OF_TRIALS)

        print("\nBEST RESULT")
        print("Value:", study.best_value)
        print("Params:", study.best_params)

        df = study.trials_dataframe()
        df.to_csv(EXPERIMENT_DIR / "optuna_results.csv", index=False, mode="w")
        # for t in study.trials:
        #     print(
        #         f"Trial {t.number:2d} | "
        #         f"value={t.value:.6f} | "
        #         f"params={t.params} | "
        #         f"state={t.state}"
        #     )

    finally:
        # Restore the original content of the bdm.toml and params.json files always, even if an error occurs during the optimization process
        BDM_TOML_PATH.write_text(original_bdm)
        PARAMS_PATH.write_text(original_params)
