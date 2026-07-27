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
EXPERIMENT_ID = 16
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