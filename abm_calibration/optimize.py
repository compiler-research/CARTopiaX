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

The simulations run inside this Python process through the cppyy bundled with
BioDynaMo's ROOT, and their output is read in memory. Build the simulation
first, then source BioDynaMo once and run the script:
    source <path-to-biodynamo>/bin/thisbdm.sh
    python3 abm_calibration/optimize.py
"""

import ctypes
import json
import logging
import os
from pathlib import Path

import numpy as np
import optuna
import pandas as pd

# Change this: File Parameters for the desired experiment
EXPERIMENT_ID = 0
# "total": compares against target_data/final_data.csv
# "radius": compares against target_data/data_dependent_on_radius_tumor.csv
# "custom": implement your own error function in compute_custom_error()
MODE = "total"

SEED = 42
# You can set the number of trials to 0 to skip the optimization and just load the best result from the database
NUMBER_OF_TRIALS = 8
# Number of Monte Carlo simulations to run for each trial. Use one for an aproximation of the error with a single montecarlo run
NUMBER_MONTE_CARLO = 3
# Minutes simulated in each run
TOTAL_MINUTES_TO_SIMULATE = 1440
# Minute at which the radial profile is compared in "radius" and "custom" modes (4320 = 72 hours).
# It must be an output minute of the simulation (every 12 hours by default) and not later than TOTAL_MINUTES_TO_SIMULATE
RADIUS_PROFILE_MINUTE = 4320

# Other Hyperparameters for the optimization process
# LOGGING
logging.basicConfig(level=logging.INFO)

# Root directory of the repository
PROJECT_ROOT = Path(__file__).resolve().parent.parent
# Simulation library built with the project
LIB_PATH = PROJECT_ROOT / "build" / "libCARTopiaX.so"

# Directory with the target data
TARGET_DIR = Path(__file__).resolve().parent / "target_data"

# Directory for this experiment
EXPERIMENT_DIR = PROJECT_ROOT / "abm_calibration" / f"experiment_{EXPERIMENT_ID}"
EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)

#######################################################
# Load the simulation in this process (only once)
#######################################################

BDMSYS = os.environ.get("BDMSYS")
if not BDMSYS:
    raise RuntimeError(
        "BioDynaMo is not sourced. Run once in your shell before this script:\n"
        "  source <path-to-biodynamo>/bin/thisbdm.sh"
    )
if not LIB_PATH.exists():
    raise FileNotFoundError(
        f"Simulation library not found at {LIB_PATH}. Build the simulation first."
    )


ctypes.CDLL("libstdc++.so.6", mode=ctypes.RTLD_GLOBAL)
import cppyy  # noqa: E402
for lib in ("libbiodynamo.so", "libbiodynamo-pv.so"):
    lib_path = Path(BDMSYS) / "lib" / lib
    if lib_path.exists():
        ctypes.CDLL(str(lib_path), mode=ctypes.RTLD_GLOBAL)

for include_dir in (Path(BDMSYS) / "include", PROJECT_ROOT / "src"):
    if include_dir.is_dir():
        cppyy.add_include_path(str(include_dir))

cppyy.load_library(str(LIB_PATH))
cppyy.include("cart_tumor.h")
bdm = cppyy.gbl.bdm
os.chdir(PROJECT_ROOT)


#######################################################
# Change this: Simulations to run the ABM with the given parameters
#######################################################


# Function to run the ABM with the given parameters. Returns its output in memory
def run_ABM(params, seed):
    # Change this: parameter to be optimized in the ABM simulation
    initial_oxygen_level = params["initial_oxygen_level"]
    default_oxygen_consumption_tumor_cell = params[
        "default_oxygen_consumption_tumor_cell"
    ]

    # Change this configuration for the ABM run (same keys as params.json)
    config = {
        "seed": seed,
        "output_performance_statistics": False,
        "total_minutes_to_simulate": TOTAL_MINUTES_TO_SIMULATE,
        "treatment": {"0": 50},
        # The data dependent on the radius is only needed by the "radius" and "custom" errors
        "output_information_dependent_on_radius": MODE in ("radius", "custom"),
        "initial_oxygen_level": initial_oxygen_level,
        "default_oxygen_consumption_tumor_cell": default_oxygen_consumption_tumor_cell,
    }

    # Load the parameters as the simulation does with params.json, so the values
    # derived from them are computed in the same way
    sim_param = bdm.SimParam()
    sim_param.LoadParamsFromString(json.dumps(config))

    # Run the ABM simulation in this process (without CSV files nor visualization),
    # its output is stored in sim_result
    sim_result = bdm.SimulationResult()
    if bdm.RunSimulation(sim_param, sim_result) != 0:
        raise RuntimeError(f"Simulation failed with config {config}")
    return sim_result


#######################################################
# Change this: Error functions (examples)
#######################################################


# Error using the data of the whole tumor as a function of time (final_data.csv)
def compute_error_total(sim_result):
    metric = "average_oxygen_all_cells"  # Change this: choose the metric to compare, e.g., "average_oxygen_all_cells", "tumor_radius", etc.
    return mse_against_target(sim_result, "final_data.csv", [metric])


# Error using the data as a function of the radius at a fixed time (data_dependent_on_radius_tumor.csv)
def compute_error_radius(sim_result):
    csv_name = "data_dependent_on_radius_tumor.csv"
    metric = "tumor_cells_type5_dead"  # Change this: prefix of the metric columns, e.g., "tumor_cells_type5_dead"

    df_t = read_target(csv_name)
    if df_t is None:
        return float("inf")

    columns = [c for c in df_t.columns if c.startswith(f"{metric}_radius_")]
    if not columns:
        logging.error("No columns for metric %s in %s", metric, TARGET_DIR / csv_name)
        return float("inf")
    return mse_against_target(
        sim_result, csv_name, columns, minute=RADIUS_PROFILE_MINUTE
    )


# Change this: write your own error. Example: dead cells at the center of a tumor
# of 3000 micrometers of radius split into 20 measurement rings at RADIUS_PROFILE_MINUTE
def compute_custom_error(sim_result):
    column = "tumor_cells_type5_dead_radius_0_to_150"
    df_s = simulation_output_to_df(sim_result.radius_data, ["total_minutes", column])
    if df_s is None:
        return float("inf")

    rows = df_s[df_s["total_minutes"] == RADIUS_PROFILE_MINUTE]
    if rows.empty:
        logging.error("No simulation output at minute %s", RADIUS_PROFILE_MINUTE)
        return float("inf")
    target_value_center = 21
    return abs(float(rows[column].iloc[0]) - target_value_center)


#######################################################
# Change this: Objective function for the Optuna optimization process
#######################################################


def objective(trial):
    # Change this: Define the parameters to be optimized and their ranges
    params = {
        "initial_oxygen_level": trial.suggest_float("initial_oxygen_level", 30, 40),
        "default_oxygen_consumption_tumor_cell": trial.suggest_float(
            "default_oxygen_consumption_tumor_cell", 7, 14
        ),
    }

    logging.info(f"Trial {trial.number} | params={params}")

    # Compute the error as the average of the errors from multiple Monte Carlo simulations varying the seed
    total_error = 0
    for seed in np.random.randint(0, 10000, NUMBER_MONTE_CARLO):
        sim_result = run_ABM(params, int(seed))

        error = compute_error(sim_result)
        total_error += error

    error = total_error / NUMBER_MONTE_CARLO

    logging.info(f"Trial {trial.number} | error={error}")

    return error


#######################################################
# Auxiliary functions
#######################################################


# Function to compute the error between the ABM simulation results and the target data
def compute_error(sim_result):
    if MODE == "total":
        return compute_error_total(sim_result)
    if MODE == "radius":
        return compute_error_radius(sim_result)
    if MODE == "custom":
        return compute_custom_error(sim_result)
    raise ValueError(f"Unknown MODE: {MODE}")


# Target CSVs, read only once because they are the same for every run
_TARGET_CACHE = {}


# Function to read a target CSV (None if it does not exist)
def read_target(csv_name):
    if csv_name not in _TARGET_CACHE:
        target = TARGET_DIR / csv_name
        if not target.exists():
            logging.error("Missing CSV: %s", target)
            return None
        _TARGET_CACHE[csv_name] = pd.read_csv(target)
    return _TARGET_CACHE[csv_name]


# Function to convert the in-memory output of the simulation (column name -> values)
# into a DataFrame with the given columns (None if a column is missing)
def simulation_output_to_df(output, columns):
    data = {}
    for column in columns:
        # count() instead of [] because std::map::operator[] inserts missing keys
        if output.count(column) == 0:
            logging.error("Missing simulation column: %s", column)
            return None
        values = output[column]
        data[column] = np.fromiter(values, dtype=np.float64, count=len(values))
    return pd.DataFrame(data)


# MSE between simulation and target for the given csv file and columns, merged on total_minutes
# (if minute is given, only that time is compared)
def mse_against_target(sim_result, csv_name, columns, minute=None):
    df_t = read_target(csv_name)
    if df_t is None:
        return float("inf")
    df_t = df_t[["total_minutes", *columns]]

    # Simulation output with the same columns as the CSV file written in output/
    output = (
        sim_result.final_data
        if csv_name == "final_data.csv"
        else sim_result.radius_data
    )
    df_s = simulation_output_to_df(output, ["total_minutes", *columns])
    if df_s is None:
        return float("inf")

    if minute is not None:
        df_t = df_t[df_t["total_minutes"] == minute]
        df_s = df_s[df_s["total_minutes"] == minute]

    merged = pd.merge(
        df_t, df_s, on="total_minutes", how="inner", suffixes=("_t", "_s")
    )
    if merged.empty:
        logging.error("No common minutes between target and simulation")
        return float("inf")

    squared = [(merged[f"{c}_s"] - merged[f"{c}_t"]) ** 2 for c in columns]
    return float(pd.concat(squared, axis=1).stack().mean())


#######################################################
# Main
#######################################################

if __name__ == "__main__":
    # The radial profile can only be compared at a minute that is simulated
    if MODE in ("radius", "custom") and RADIUS_PROFILE_MINUTE > TOTAL_MINUTES_TO_SIMULATE:
        raise ValueError(
            f"RADIUS_PROFILE_MINUTE ({RADIUS_PROFILE_MINUTE}) is later than "
            f"TOTAL_MINUTES_TO_SIMULATE ({TOTAL_MINUTES_TO_SIMULATE})"
        )

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

    # A crash inside a simulation stops the whole process and leaves its trial as
    # RUNNING. This script runs one trial at a time, so trials still RUNNING come
    # from a previous execution that was stopped and are marked as failed
    for stale_trial in study.get_trials(
        deepcopy=False, states=(optuna.trial.TrialState.RUNNING,)
    ):
        study.tell(stale_trial.number, state=optuna.trial.TrialState.FAIL)

    # The visualization export is disabled for each run by RunSimulation, so
    # bdm.toml does not need to be modified
    study.optimize(objective, n_trials=NUMBER_OF_TRIALS)

    print("\nBEST RESULT")
    print("Value:", study.best_value)
    print("Params:", study.best_params)

    df = study.trials_dataframe()
    df.to_csv(EXPERIMENT_DIR / "optuna_results.csv", index=False, mode="w")
