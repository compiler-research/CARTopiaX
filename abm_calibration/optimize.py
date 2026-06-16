import json
import optuna
import numpy as np
import logging
from pathlib import Path
import os
import subprocess

# Change this: File Parameters for the desired experiment
EXPERIMENT_ID =  1
SEED = 42
#You can set the number of trials to 0 to skip the optimization and just load the best result from the database
NUMBER_OF_TRIALS = 1
# Number of Monte Carlo simulations to run for each trial. Use one for an aproximation of the error with a single montecarlo run
NUMBER_MONTE_CARLO = 3
#BioDynaMo directory to execute the comand source thisbdm.sh, you can change it to your own path
BIODYNAMO_DIR = " /home/usuario/Desktop/biodynamo/build/bin/thisbdm.sh"

# Other Hyperparameters for the optimization process
# LOGGING
logging.basicConfig(level=logging.INFO)

# Directory for the visualization configuration file bdm.toml
BDM_TOML_PATH = Path("./bdm.toml")
# Directory for the parameters file params.json
PARAMS_PATH = Path(__file__).resolve().parent.parent / "params.json"

#Directory for this experiment
EXPERIMENT_DIR = Path("abm_calibration") / f"experiment_{EXPERIMENT_ID}"
EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)


#Function to run the ABM with the given parameters
def run_ABM(params, seed):
    # Change this: parameter to be optimized in the ABM simulation 
    kill_rate_cart = params["kill_rate_cart"]
    adhesion_rate_cart = params["adhesion_rate_cart"]

   # Change this configuration for the ABM run
    config = {
        "seed": seed,
        "output_performance_statistics": False,
        "total_minutes_to_simulate": 1440,
        "initial_tumor_radius": 40.0,
        "treatment": {
            "0": 50
        },
        "kill_rate_cart": kill_rate_cart,
        "adhesion_rate_cart": adhesion_rate_cart
    }

    # Save the config parameters for the run to the params.json file
    with open(PARAMS_PATH, "w") as f:
        json.dump(config, f, indent=2)


    # Load the ByoDynaMo environment and run the ABM simulation using the BioDynaMo executable
    subprocess.run(
        [
            "bash",
            "-c",
            f"source {BIODYNAMO_DIR} && bdm run"
        ],
        check=True
    )

# Change This: Function to compute the error between the ABM simulation results and the experimental data
def compute_error():
    return 0.0  # Placeholder for the actual error computation logic

# Objective function for the Optuna optimization process
def objective(trial):
   # Change this: Define the parameters to be optimized and their ranges 
    params = {
        "kill_rate_cart": trial.suggest_float("kill_rate_cart", 0.0, 0.1),
        "adhesion_rate_cart": trial.suggest_float("adhesion_rate_cart", 0.0, 0.1),
    }

    logging.info(f"Trial {trial.number} | params={params}")

    # Compute the error as the average of the errors from multiple Monte Carlo simulations varying the seed
    total_error = 0
    for seed in np.random.randint(0, 10000, NUMBER_MONTE_CARLO):
        run_ABM(params, int(seed))
        
        error = compute_error()
        total_error += error

    error= total_error / NUMBER_MONTE_CARLO

    logging.info(f"Trial {trial.number} | error={error}")

    return error

# Fix the random seed for reproducibility
np.random.seed(SEED)

# Create an Optuna study to optimize the parameters of the ABM
study = optuna.create_study(
    study_name="abm_calibration",
    storage = f"sqlite:///{EXPERIMENT_DIR / 'abm_optuna.db'}",
    load_if_exists=True,
    direction="minimize",
    sampler=optuna.samplers.TPESampler(seed=SEED)
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

        study.optimize(
            objective,
            n_trials=NUMBER_OF_TRIALS
        )

        print("\nBEST RESULT")
        print("Value:", study.best_value)
        print("Params:", study.best_params)

        df = study.trials_dataframe()
        df.to_csv(EXPERIMENT_DIR / "optuna_results.csv", index=False, mode="w")

    finally:
        # Restore the original content of the bdm.toml and params.json files always, even if an error occurs during the optimization process
        BDM_TOML_PATH.write_text(original_bdm)
        PARAMS_PATH.write_text(original_params)
