import json
import logging
import subprocess
from pathlib import Path

import numpy as np
import optuna
import pandas as pd

# Change this: File Parameters for the desired experiment
EXPERIMENT_ID = 1
SEED = 42
# You can set the number of trials to 0 to skip the optimization and just load the best result from the database
NUMBER_OF_TRIALS = 8
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
    initial_oxygen_level = params["initial_oxygen_level"]
    default_oxygen_consumption_tumor_cell = params["default_oxygen_consumption_tumor_cell"]

    # Change this configuration for the ABM run
    config = {
        "seed": seed,
        "output_performance_statistics": False,
        "total_minutes_to_simulate": 1440,
        "initial_tumor_radius": 40.0,
        "treatment": {"0": 50},
        "initial_oxygen_level": initial_oxygen_level,
        "default_oxygen_consumption_tumor_cell": default_oxygen_consumption_tumor_cell,
    }

    # Save the config parameters for the run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)

    # Load the ByoDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(["bash", "-c", f"source {BIODYNAMO_DIR} && bdm run"], check=True)


# Change This: Function to compute the error between the ABM simulation results and the experimental data
def compute_error():
    target = Path(__file__).resolve().parent / "target_data" / "final_data.csv"
    sim_dir = Path(__file__).resolve().parent.parent / "output" / "final_data.csv"

    if not target.exists():
        logging.error("Missing target CSV: %s", target)
        return float("inf")

    if not sim_dir.exists():
        logging.error("Missing simulation CSV: %s", sim_dir)
        return float("inf")

    df_t = pd.read_csv(target, usecols=["total_minutes", "average_oxygen_cancer_cells"])
    df_s = pd.read_csv(sim_dir, usecols=["total_minutes", "average_oxygen_cancer_cells"])

    # Merge the target and simulation dataframes on the "total_minutes" column to align the data for error computation
    merged = pd.merge(
        df_t, df_s, on="total_minutes", how="inner", suffixes=("_t", "_s")
    )

    if merged.empty:
        logging.error("No common minutes between target and simulation")
        return float("inf")

    # Compute the mean squared error (MSE) between the average oxygen levels in the target and simulation data
    mse = ((merged["average_oxygen_cancer_cells_s"] - merged["average_oxygen_cancer_cells_t"]) ** 2).mean()
    return float(mse)


# Objective function for the Optuna optimization process
def objective(trial):
    # Change this: Define the parameters to be optimized and their ranges
    params = {
        "initial_oxygen_level": trial.suggest_float("initial_oxygen_level", 30, 40),
        "default_oxygen_consumption_tumor_cell": trial.suggest_float("default_oxygen_consumption_tumor_cell", 7, 14),
    }

    logging.info(f"Trial {trial.number} | params={params}")

    # Compute the error as the average of the errors from multiple Monte Carlo simulations varying the seed
    total_error = 0
    for seed in np.random.randint(0, 10000, NUMBER_MONTE_CARLO):
        run_ABM(params, int(seed))

        error = compute_error()
        total_error += error

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

    finally:
        # Restore the original content of the bdm.toml and params.json files always, even if an error occurs during the optimization process
        BDM_TOML_PATH.write_text(original_bdm)
        PARAMS_PATH.write_text(original_params)
