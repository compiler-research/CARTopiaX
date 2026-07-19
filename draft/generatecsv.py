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
# BioDynaMo directory to execute the comand source thisbdm.sh, you can change it to your own path
BIODYNAMO_DIR = " /home/usuario/Desktop/biodynamo/build/bin/thisbdm.sh"
EXPERIMENT_DIR = Path("abm_calibration") / f"experiment_{EXPERIMENT_ID}"

# Create an Optuna study to optimize the parameters of the ABM
study = optuna.create_study(
    study_name="abm_calibration",
    storage=f"sqlite:///{EXPERIMENT_DIR / 'abm_optuna.db'}",
    load_if_exists=True,
    direction="minimize",
    sampler=optuna.samplers.TPESampler(seed=SEED),
)

df = study.trials_dataframe()
df.to_csv(EXPERIMENT_DIR / "optuna_results.csv", index=False, mode="w")