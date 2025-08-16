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
import sys

# Configuration flags
SAVE_TO_FILE = False  # Set to True to save output to file, False to print to console
OUTPUT_FILENAME = "comparison_results_freeOxygen_freeOncoprotein_freeRandomRate.txt"  # Name of the output file
num_files = 60  # Num
# Redirect output if flag is True
if SAVE_TO_FILE:
    original_stdout = sys.stdout
    sys.stdout = open(OUTPUT_FILENAME, 'w')

for i in range(num_files):
    # Read data from new CSV format with _processed suffix
    df_cartopiax = pd.read_csv('processed/simulation_data_CARTopiaX'+str(i)+'_processed.csv')
    df_paper = pd.read_csv('processed/simulation_data'+str(i)+'_processed.csv')
    
    # Extract values for CARTopiaX data
    min_cartopiax = {
        'oxygen_level': df_cartopiax['min_oxygen_level'].iloc[0],
        'oncoproteine_level': df_cartopiax['min_oncoproteine_level'].iloc[0],
        'base_transition_rate': df_cartopiax['min_base_transition_rate'].iloc[0],
        'final_rate_transition': df_cartopiax['min_final_rate_transition'].iloc[0],
        'probability_necrosis': df_cartopiax['min_probability_necrosis'].iloc[0]
    }
    
    mean_cartopiax = {
        'oxygen_level': df_cartopiax['avg_oxygen_level'].iloc[0],
        'oncoproteine_level': df_cartopiax['avg_oncoproteine_level'].iloc[0],
        'base_transition_rate': df_cartopiax['avg_base_transition_rate'].iloc[0],
        'final_rate_transition': df_cartopiax['avg_final_rate_transition'].iloc[0],
        'probability_necrosis': df_cartopiax['avg_probability_necrosis'].iloc[0]
    }
    
    max_cartopiax = {
        'oxygen_level': df_cartopiax['max_oxygen_level'].iloc[0],
        'oncoproteine_level': df_cartopiax['max_oncoproteine_level'].iloc[0],
        'base_transition_rate': df_cartopiax['max_base_transition_rate'].iloc[0],
        'final_rate_transition': df_cartopiax['max_final_rate_transition'].iloc[0],
        'probability_necrosis': df_cartopiax['max_probability_necrosis'].iloc[0]
    }
    
    # Extract values for the their model's data
    min_paper = {
        'oxygen_level': df_paper['min_oxygen_level'].iloc[0],
        'oncoprotein_level': df_paper['min_oncoproteine_level'].iloc[0],
        'base_transition_rate': df_paper['min_base_transition_rate'].iloc[0],
        'final_rate_transition': df_paper['min_final_rate_transition'].iloc[0],
        'probability_necrosis': df_paper['min_probability_necrosis'].iloc[0]
    }
    
    mean_paper = {
        'oxygen_level': df_paper['avg_oxygen_level'].iloc[0],
        'oncoprotein_level': df_paper['avg_oncoproteine_level'].iloc[0],
        'base_transition_rate': df_paper['avg_base_transition_rate'].iloc[0],
        'final_rate_transition': df_paper['avg_final_rate_transition'].iloc[0],
        'probability_necrosis': df_paper['avg_probability_necrosis'].iloc[0]
    }
    
    max_paper = {
        'oxygen_level': df_paper['max_oxygen_level'].iloc[0],
        'oncoprotein_level': df_paper['max_oncoproteine_level'].iloc[0],
        'base_transition_rate': df_paper['max_base_transition_rate'].iloc[0],
        'final_rate_transition': df_paper['max_final_rate_transition'].iloc[0],
        'probability_necrosis': df_paper['max_probability_necrosis'].iloc[0]
    }

    # Print as table
    print(f"\nComparison for simulation {i}:")
    print(f"{'Metric':<25} {'CARTopiaX Min':>12} {'CARTopiaX Avg':>12} {'CARTopiaX Max':>12} {'Paper Min':>12} {'Paper Avg':>12} {'Paper Max':>12}")
    print("-" * 97)
    print(f"{'Oxygen level':<25} {min_cartopiax['oxygen_level']:>12.6f} {mean_cartopiax['oxygen_level']:>12.6f} {max_cartopiax['oxygen_level']:>12.6f} {min_paper['oxygen_level']:>12.6f} {mean_paper['oxygen_level']:>12.6f} {max_paper['oxygen_level']:>12.6f}")
    print(f"{'Oncoprotein level':<25} {min_cartopiax['oncoprotein_level']:>12.6f} {mean_cartopiax['oncoprotein_level']:>12.6f} {max_cartopiax['oncoprotein_level']:>12.6f} {min_paper['oncoprotein_level']:>12.6f} {mean_paper['oncoprotein_level']:>12.6f} {max_paper['oncoprotein_level']:>12.6f}")
    print(f"{'Base transition rate':<25} {min_cartopiax['base_transition_rate']:>12.6f} {mean_cartopiax['base_transition_rate']:>12.6f} {max_cartopiax['base_transition_rate']:>12.6f} {min_paper['base_transition_rate']:>12.6f} {mean_paper['base_transition_rate']:>12.6f} {max_paper['base_transition_rate']:>12.6f}")
    print(f"{'Transition rate':<25} {min_cartopiax['final_rate_transition']:>12.6f} {mean_cartopiax['final_rate_transition']:>12.6f} {max_cartopiax['final_rate_transition']:>12.6f} {min_paper['final_rate_transition']:>12.6f} {mean_paper['final_rate_transition']:>12.6f} {max_paper['final_rate_transition']:>12.6f}")
    print(f"{'Necrosis probability':<25} {min_cartopiax['probability_necrosis']:>12.6f} {mean_cartopiax['probability_necrosis']:>12.6f} {max_cartopiax['probability_necrosis']:>12.6f} {min_paper['probability_necrosis']:>12.6f} {mean_paper['probability_necrosis']:>12.6f} {max_paper['probability_necrosis']:>12.6f}")

# Restore original stdout if output was redirected to file
if SAVE_TO_FILE:
    sys.stdout.close()
    sys.stdout = original_stdout
    print(f"Output saved to {OUTPUT_FILENAME}")
