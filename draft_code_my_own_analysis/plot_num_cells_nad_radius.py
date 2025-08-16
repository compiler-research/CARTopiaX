"""
Copyright 2025 compiler-research.org, Salvador de la Torre Gonzalez, Luciana Melina Luque

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

This file contains a model developed under Google Summer of Code (GSoC)
for the compiler-research.org organization.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Paths to the data files
csv_data_path = 'final_data.csv'
dat_data_path = 'DatosFinales.dat'
days = 30.1

# Read the CSV data
df_csv = pd.read_csv(csv_data_path)

column_names = ['#tiempo', 'volumen', 'volumen2', 'radio', 'celulas tumorales', 
                'celulas muertas', 'tumor_muerto', 'todas las celulas']
df_dat = pd.read_csv(dat_data_path, sep=' ', names=column_names, skiprows=1)

df_dat['#tiempo'] = pd.to_numeric(df_dat['#tiempo'])

df_dat['total_days'] = df_dat['#tiempo'] / 1440  # 1440 minutes in a day

filtered_df_csv = df_csv[df_csv['total_days'] <= days]
filtered_df_dat = df_dat[df_dat['total_days'] <= days]

fig, ax1 = plt.subplots(figsize=(10, 6))

color1 = 'blue'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Tumor Radius', color=color1)
ax1.plot(filtered_df_csv['total_days'], filtered_df_csv['tumor_radius'], 
         color=color1, marker='o', label='Tumor Radius: CARTopiaX')
ax1.plot(filtered_df_dat['total_days'], filtered_df_dat['radio'], 
         color=color1, linestyle='--', marker='x', label='Tumor Radius: Model from Nature paper')
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()
color2 = 'red'
ax2.set_ylabel('Number of Cells', color=color2)
ax2.plot(filtered_df_csv['total_days'], filtered_df_csv['num_cells'], 
         color=color2, marker='s', label='Number of Cells: Model in BioDynaMo')
ax2.plot(filtered_df_dat['total_days'], filtered_df_dat['todas las celulas'], 
         color=color2, linestyle='--', marker='+', label='Number of Cells: Model from Nature paper')
ax2.tick_params(axis='y', labelcolor=color2)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

plt.title('Tumor Radius and Cell Count over Time')

plt.tight_layout()
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('./no_cart_free_oxygen_free_oncoproteine_free_random_rate_tumor_radius_and_cells_comparison.png', dpi=300)
print("Plot saved")