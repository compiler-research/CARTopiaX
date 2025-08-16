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
import numpy as np

# Configuration
num_files = 60  # Files from 0 to 59
time_interval = 0.5  # Each file represents half a day

# Initialize lists to store data
days = []
oxygen_mine = []
oxygen_paper = []
oncoproteine_mine = []
oncoproteine_paper = []

# Read data from all files
for i in range(num_files):
    # Calculate the day for this file
    day = i * time_interval
    days.append(day)
    
    try:
        # Read data from processed CSV files
        df_mine = pd.read_csv(f'processed/simulation_data_CARTopiaX{i}_processed.csv')
        df_paper = pd.read_csv(f'processed/simulation_data{i}_processed.csv')
        
        # Extract average values for plotting
        oxygen_mine.append(df_mine['avg_oxygen_level'].iloc[0])
        oxygen_paper.append(df_paper['avg_oxygen_level'].iloc[0])
        oncoproteine_mine.append(df_mine['avg_oncoproteine_level'].iloc[0])
        oncoproteine_paper.append(df_paper['avg_oncoproteine_level'].iloc[0])
        
    except FileNotFoundError:
        print(f"Warning: File for simulation {i} not found, skipping...")
        # Remove the day we just added if file doesn't exist
        days.pop()

# Create the plot
fig, ax1 = plt.subplots(figsize=(12, 8))

# Plot oxygen levels on the primary y-axis
color1 = 'blue'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Average Oxygen Level', color=color1)
ax1.plot(days, oxygen_mine, color=color1, marker='o', markersize=4, 
         label='Average Oxygen Level: Model in BioDynaMo', linewidth=2)
ax1.plot(days, oxygen_paper, color=color1, linestyle='--', marker='x', markersize=4,
         label='Average Oxygen Level: Model from Nature paper', linewidth=2)
ax1.tick_params(axis='y', labelcolor=color1)

# Create secondary y-axis for oncoproteine levels
ax2 = ax1.twinx()
color2 = 'red'
ax2.set_ylabel('Average Oncoprotein Level', color=color2)
ax2.plot(days, oncoproteine_mine, color=color2, marker='s', markersize=4,
         label='Average Oncoprotein Level: Model in BioDynaMo', linewidth=2)
ax2.plot(days, oncoproteine_paper, color=color2, linestyle='--', marker='+', markersize=4,
         label='Average Oncoprotein Level: Model from Nature paper', linewidth=2)
ax2.tick_params(axis='y', labelcolor=color2)

# Combine legends from both axes and position at top center
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2)

# Set title and formatting
plt.title('Oxygen and Oncoprotein Levels over Time Comparison', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.grid(True, linestyle='--', alpha=0.7)

# Save the plot
plt.savefig('./no_cart_oxygen_oncoproteine_comparison.png', dpi=300, bbox_inches='tight')
print("Plot saved as 'no_cart_oxygen_oncoproteine_comparison.png'")
plt.show()

# Print summary statistics
print(f"\nSummary Statistics:")
print(f"Time range: {min(days):.1f} to {max(days):.1f} days")
print(f"Oxygen - Model in BioDynaMo: min={min(oxygen_mine):.6f}, max={max(oxygen_mine):.6f}, avg={np.mean(oxygen_mine):.6f}")
print(f"Oxygen - Model from Nature paper: min={min(oxygen_paper):.6f}, max={max(oxygen_paper):.6f}, avg={np.mean(oxygen_paper):.6f}")
print(f"Oncoproteine - Model in BioDynaMo: min={min(oncoproteine_mine):.6f}, max={max(oncoproteine_mine):.6f}, avg={np.mean(oncoproteine_mine):.6f}")
print(f"Oncoproteine - Model from Nature paper: min={min(oncoproteine_paper):.6f}, max={max(oncoproteine_paper):.6f}, avg={np.mean(oncoproteine_paper):.6f}")