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
import matplotlib.pyplot as plt

# Read CSV files
df_mine = pd.read_csv('forces.csv', header=None)
df_mine_positive = df_mine[df_mine.iloc[:, 1] >= 0]

df_theirs = pd.read_csv('fuerzas.csv', header=None)
df_theirs.iloc[:, 0] = df_theirs.iloc[:, 0] - 0.1
df_theirs_positive = df_theirs[df_theirs.iloc[:, 1] >= 0]
# Define column names for clarity
# Column 0: time, Columns 1-3: coordinate_x, coordinate_y, coordinate_z
# Columns 4-6: not relevant, Column 7: displacement_norm

# Create figure with 4 subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Time vs Coordinate X
ax1.plot(df_mine_positive.iloc[:, 0], df_mine_positive.iloc[:, 1], 
         color='blue', label='My data', linewidth=2, linestyle='-')
ax1.plot(df_theirs_positive.iloc[:, 0], df_theirs_positive.iloc[:, 1], 
         color='red', label='Their data', linewidth=2, linestyle='--')
ax1.set_xlabel('Time')
ax1.set_ylabel('Coordinate X')
ax1.set_title('Time vs Coordinate X')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Time vs Coordinate Y
ax2.plot(df_mine_positive.iloc[:, 0], df_mine_positive.iloc[:, 2], 
         color='blue', label='My data', linewidth=2, linestyle='-')
ax2.plot(df_theirs_positive.iloc[:, 0], df_theirs_positive.iloc[:, 2], 
         color='red', label='Their data', linewidth=2, linestyle='--')
ax2.set_xlabel('Time')
ax2.set_ylabel('Coordinate Y')
ax2.set_title('Time vs Coordinate Y')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Time vs Coordinate Z
ax3.plot(df_mine_positive.iloc[:, 0], df_mine_positive.iloc[:, 3], 
         color='blue', label='My data', linewidth=2, linestyle='-')
ax3.plot(df_theirs_positive.iloc[:, 0], df_theirs_positive.iloc[:, 3], 
         color='red', label='Their data', linewidth=2, linestyle='--')
ax3.set_xlabel('Time')
ax3.set_ylabel('Coordinate Z')
ax3.set_title('Time vs Coordinate Z')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Time vs Displacement Norm
ax4.plot(df_mine_positive.iloc[:, 0], df_mine_positive.iloc[:, 7], 
         color='blue', label='My data', linewidth=2, linestyle='-')
ax4.plot(df_theirs_positive.iloc[:, 0], df_theirs_positive.iloc[:, 7], 
         color='red', label='Their data', linewidth=2, linestyle='--')
ax4.set_xlabel('Time')
ax4.set_ylabel('Displacement Norm')
ax4.set_title('Time vs Displacement Norm')
ax4.legend()
ax4.grid(True, alpha=0.3)

# Adjust layout to prevent overlapping
plt.tight_layout()

# Save the figure
plt.savefig('force_comparison_plots_fixed_cell_volume_adhesion.png', dpi=300, bbox_inches='tight')
print("Plot saved")

