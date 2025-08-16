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
import os

name = 'simulation_data' #simulation_data_mine #simulation_data_theirs
num_files = 60
def process_csv_file(input_filename, output_filename):
    """
    Process a CSV file to calculate min, avg, max for each column
    and save to a new CSV with '_processed' suffix
    """
    # Read the original CSV
    df = pd.read_csv(input_filename, 
                    names=['oxygen_level', 'oncoproteine_level', 'base_transition_rate',
                           'final_rate_transition', 'probability_necrosis'])
    
    # Calculate statistics
    min_vals = df.min()
    mean_vals = df.mean()
    max_vals = df.max()
    
    # Create processed data with 15 columns (min, avg, max for each of 5 metrics)
    processed_data = {
        'min_oxygen_level': [min_vals['oxygen_level']],
        'avg_oxygen_level': [mean_vals['oxygen_level']],
        'max_oxygen_level': [max_vals['oxygen_level']],
        'min_oncoproteine_level': [min_vals['oncoproteine_level']],
        'avg_oncoproteine_level': [mean_vals['oncoproteine_level']],
        'max_oncoproteine_level': [max_vals['oncoproteine_level']],
        'min_base_transition_rate': [min_vals['base_transition_rate']],
        'avg_base_transition_rate': [mean_vals['base_transition_rate']],
        'max_base_transition_rate': [max_vals['base_transition_rate']],
        'min_final_rate_transition': [min_vals['final_rate_transition']],
        'avg_final_rate_transition': [mean_vals['final_rate_transition']],
        'max_final_rate_transition': [max_vals['final_rate_transition']],
        'min_probability_necrosis': [min_vals['probability_necrosis']],
        'avg_probability_necrosis': [mean_vals['probability_necrosis']],
        'max_probability_necrosis': [max_vals['probability_necrosis']]
    }
    
    # Create DataFrame and save
    processed_df = pd.DataFrame(processed_data)
    processed_df.to_csv(output_filename, index=False)
    print(f"Processed {input_filename} -> {output_filename}")

def read_processed_csv(filename):
    """
    Read a processed CSV file and return the DataFrame
    """
    return pd.read_csv(filename)

def main():
    for i in range(num_files):
        input_file = f'{name}{i}.csv'
        output_file = f'processed/{name}{i}_processed.csv'

        if os.path.exists(input_file):
            process_csv_file(input_file, output_file)
        else:
            print(f"Warning: {input_file} not found")
    
    print("\nProcessing complete!")
    
if __name__ == "__main__":
    main()