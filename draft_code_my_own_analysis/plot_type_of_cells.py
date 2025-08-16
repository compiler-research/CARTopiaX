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

# Read the CSV file
df = pd.read_csv('./final_data.csv')

# Calculate percentages for our model data
df['type1_percentage'] = 100 * df['tumor_cells_type1'] / df['num_tumor_cells']
df['type2_percentage'] = 100 * df['tumor_cells_type2'] / df['num_tumor_cells']
df['type3_percentage'] = 100 * df['tumor_cells_type3'] / df['num_tumor_cells']
df['type4_percentage'] = 100 * df['tumor_cells_type4'] / df['num_tumor_cells']

# Nature Paper data - days and data processing
days = [0,719,1439,2160,2880,3600,4320,5040,5760,6480,7200,7920,8640,9360,10080,10800,11520,12240,12960,13680,14400,15120,15840,16560,17280,18000,18719,19439,20159,20879,21599,22319,23039,23759,24479,25199,25919,26639,27359,28079,28799,29519,30239,30959,31679,32399,33119,33839,34559,35279,35999,36719,37439,38159,38879,39599,40319,41039,41759,42479,43199]

quantities_type1 = []
quantities_type2 = []
quantities_type3 = []
quantities_type4 = []
time_days = []

for i in range(0, len(days)):            
    file_path = 'Datos_{}.xyz'.format(days[i])            
    
    # Read data file
    data = pd.read_csv('./'+file_path+'', sep=" ", skiprows=1, header=None)
    data.columns = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q"]

    count_type1 = 0
    count_type2 = 0
    count_type3 = 0
    count_type4 = 0
    count_type0 = 0

    # Count cells by type
    for j in range(0, len(data)):
        if data.m[j] == 1:
            count_type1 += 1            
        elif data.m[j] == 2:
            count_type2 += 1            
        elif data.m[j] == 3:
            count_type3 += 1            
        elif data.m[j] == 4:
            count_type4 += 1            
        elif data.m[j] == 0:
            count_type0 += 1

    # Calculate percentages
    total_tumor_cells = count_type1 + count_type2 + count_type3 + count_type4
    quantities_type1.append(100 * count_type1 / total_tumor_cells)
    quantities_type2.append(100 * count_type2 / total_tumor_cells)
    quantities_type3.append(100 * count_type3 / total_tumor_cells)
    quantities_type4.append(100 * count_type4 / total_tumor_cells)
    time_days.append(days[i] / 1440)  # Convert to days

# Plot the number of cells of each type over time (types 1-4)
plt.figure(figsize=(10, 6))

# Our model data (solid lines) - now as percentages
plt.plot(df['total_days'], df['type1_percentage'], color='darkred', label='Type 1 (CARTopiaX)', linewidth=2)
plt.plot(df['total_days'], df['type2_percentage'], color='red', label='Type 2 (CARTopiaX)', linewidth=2)
plt.plot(df['total_days'], df['type3_percentage'], color='darkorange', label='Type 3 (CARTopiaX)', linewidth=2)
plt.plot(df['total_days'], df['type4_percentage'], color='gold', label='Type 4 (CARTopiaX)', linewidth=2)

# Nature Paper data (dashed lines)
plt.plot(time_days, quantities_type1, color='darkred', linestyle='--', label='Type 1 (Model from Nature paper)', linewidth=2)
plt.plot(time_days, quantities_type2, color='red', linestyle='--', label='Type 2 (Model from Nature paper)', linewidth=2)
plt.plot(time_days, quantities_type3, color='darkorange', linestyle='--', label='Type 3 (Model from Nature paper)', linewidth=2)
plt.plot(time_days, quantities_type4, color='gold', linestyle='--', label='Type 4 (Model from Nature paper)', linewidth=2)

plt.xlabel('Time (days)')
plt.ylabel('Cell types [%]')
plt.title('Tumor Cell Types Comparison: no CAR-T involved')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig('./no_cart_free_oxygen_free_oncoproteine_free_random_rate_tumor_cell_types_comparison.png', dpi=300)
print("Plot saved")