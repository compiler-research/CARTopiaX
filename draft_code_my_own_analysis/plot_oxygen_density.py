# -----------------------------------------------------------------------------
# Copyright (C) 2025 Salvador de la Torre Gonzalez
# Co-author: Luciana Melina Luque
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


df_theirs = pd.read_csv("./simulation_data_theirs11.csv", 
                names=['oxygen_level', 'oncoproteine_level', 'base_transition_rate',
                        'final_rate_transition', 'probability_necrosis'])
df_mine = pd.read_csv("./simulation_data_mine12.csv", 
                names=['oxygen_level', 'oncoproteine_level', 'base_transition_rate',
                        'final_rate_transition', 'probability_necrosis'])


plt.figure(figsize=(10, 6))
# Create normalized density plots for both datasets
sns.kdeplot(df_theirs['oxygen_level'], fill=True, alpha=0.6, label='Theirs', color='blue')
sns.kdeplot(df_mine['oxygen_level'], fill=True, alpha=0.6, label='Mine', color='red')

# Calculate means
mean_theirs = df_theirs['oxygen_level'].mean()
mean_mine = df_mine['oxygen_level'].mean()

# Add vertical lines for means
plt.axvline(mean_theirs, color='blue', linestyle='--', alpha=0.8, linewidth=2)
plt.axvline(mean_mine, color='red', linestyle='--', alpha=0.8, linewidth=2)

# Add text annotations for means
plt.text(mean_theirs, plt.ylim()[1]*0.9, f'Mean (Theirs): {mean_theirs:.2f}', 
         rotation=90, verticalalignment='top', color='black', fontsize=10, fontweight='bold')
plt.text(mean_mine, plt.ylim()[1]*0.5, f'Mean (Mine): {mean_mine:.2f}', 
         rotation=90, verticalalignment='top', color='black', fontsize=10, fontweight='bold')

plt.xlabel('Oxygen Level')
plt.ylabel('Probability Density')
plt.title('Comparison of Oxygen Level')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the plot
# plt.savefig('oxygen_no_force_without_consumption.png', dpi=300, bbox_inches='tight')
plt.show()

# Compare statistics for both datasets
print("=== THEIRS ===")
print(f"Data range: {df_theirs['oxygen_level'].min():.2f} to {df_theirs['oxygen_level'].max():.2f}")
print(f"Mean: {df_theirs['oxygen_level'].mean():.2f}")
print(f"Standard deviation: {df_theirs['oxygen_level'].std():.2f}")

print("\n=== MINE ===")
print(f"Data range: {df_mine['oxygen_level'].min():.2f} to {df_mine['oxygen_level'].max():.2f}")
print(f"Mean: {df_mine['oxygen_level'].mean():.2f}")
print(f"Standard deviation: {df_mine['oxygen_level'].std():.2f}")