/*
 * Copyright 2025 compiler-research.org, Salvador de la Torre Gonzalez, Luciana
 * Melina Luque
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *     SPDX-License-Identifier: Apache-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * This file contains a model developed under Google Summer of Code (GSoC)
 * for the compiler-research.org organization.
 */

#include "hyperparams.h"
#include "core/param/param_group.h"
#include "core/real_t.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

namespace bdm {

// Define the static member kUid for SimParam
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

// Function to read a JSON file and load the parameters
// The parameters that are not found in the file keep their default
void SimParam::LoadParams(const std::string& filename) {
  nlohmann::json jfile;
  std::ifstream file(filename);
  // If the json can be opened, try to read it, otherwise leave jfile empty
  if (file.is_open()) {
    try {
      file >> jfile;
    } catch (const std::exception& e) {
      std::cerr << "Error reading JSON: " << e.what() << std::endl;
    }
  }

  // Helper functions to load parameters from JSON
  auto load_double = [&](const std::string& key, double& var) {
    if (jfile.contains(key)) {
      var = jfile[key].get<double>();
    }
  };

  auto load_int = [&](const std::string& key, int& var) {
    if (jfile.contains(key)) {
      var = jfile[key].get<int>();
    }
  };

  auto load_bool = [&](const std::string& key, bool& var) {
    if (jfile.contains(key)) {
      var = jfile[key].get<bool>();
    }
  };

  // Load parameters from JSON file
  load_int("seed", seed);
  load_bool("output_performance_statistics", output_performance_statistics);
  load_int("total_minutes_to_simulate", total_minutes_to_simulate);
  load_double("initial_tumor_radius", initial_tumor_radius);
  load_int("bounded_space_length", bounded_space_length);

  if (jfile.contains("treatment")) {
    treatment.clear();
    // Parse the JSON object to handle string keys
    for (auto& [key, value] : jfile["treatment"].items()) {
      int const day = std::stoi(key);
      int const dose = value.get<int>();
      treatment[day] = dose;
    }
  }

  load_double("dt_substances", dt_substances);
  load_double("dt_mechanics", dt_mechanics);
  load_double("dt_cycle", dt_cycle);

  if (jfile.contains("dt_step")) {
    dt_step = jfile["dt_step"].get<double>();
  } else {
    dt_step = dt_mechanics;
  }

  if (jfile.contains("output_csv_interval")) {
    output_csv_interval = jfile["output_csv_interval"].get<int>();
  } else {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    output_csv_interval = static_cast<int>(12 * 60 / dt_step);
  }

  load_double("volume_relaxation_rate_cytoplasm_apoptotic_cells",
              volume_relaxation_rate_cytoplasm_apoptotic_cells);
  load_double("volume_relaxation_rate_nucleus_apoptotic_cells",
              volume_relaxation_rate_nucleus_apoptotic_cells);
  load_double("volume_relaxation_rate_fluid_apoptotic_cells",
              volume_relaxation_rate_fluid_apoptotic_cells);

  load_double("time_apoptosis", time_apoptosis);
  load_double("reduction_consumption_dead_cells",
              reduction_consumption_dead_cells);
  load_int("resolution_grid_substances", resolution_grid_substances);

  load_double("diffusion_coefficient_oxygen", diffusion_coefficient_oxygen);
  load_double("decay_constant_oxygen", decay_constant_oxygen);
  load_double("diffusion_coefficient_immunostimulatory_factor",
              diffusion_coefficient_immunostimulatory_factor);
  load_double("decay_constant_immunostimulatory_factor",
              decay_constant_immunostimulatory_factor);
  load_double("oxygen_reference_level", oxygen_reference_level);
  load_double("initial_oxygen_level", initial_oxygen_level);
  load_double("oxygen_saturation", oxygen_saturation);

  load_double("cell_repulsion_between_tumor_tumor",
              cell_repulsion_between_tumor_tumor);
  load_double("cell_repulsion_between_cart_cart",
              cell_repulsion_between_cart_cart);
  load_double("cell_repulsion_between_cart_tumor",
              cell_repulsion_between_cart_tumor);
  load_double("cell_repulsion_between_tumor_cart",
              cell_repulsion_between_tumor_cart);

  load_double("max_relative_adhesion_distance", max_relative_adhesion_distance);
  load_double("cell_adhesion_between_tumor_tumor",
              cell_adhesion_between_tumor_tumor);
  load_double("cell_adhesion_between_cart_cart",
              cell_adhesion_between_cart_cart);
  load_double("cell_adhesion_between_cart_tumor",
              cell_adhesion_between_cart_tumor);
  load_double("cell_adhesion_between_tumor_cart",
              cell_adhesion_between_tumor_cart);
  load_int("length_box_mechanics", length_box_mechanics);

  if (jfile.contains("dnew")) {
    dnew = jfile["dnew"].get<double>();
  } else {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    dnew = 1.5 * dt_mechanics;
  }

  if (jfile.contains("dold")) {
    dold = jfile["dold"].get<double>();
  } else {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    dold = -0.5 * dt_mechanics;
  }

  load_double("rate_secretion_immunostimulatory_factor",
              rate_secretion_immunostimulatory_factor);
  load_double("saturation_density_immunostimulatory_factor",
              saturation_density_immunostimulatory_factor);
  load_double("oncoprotein_mean", oncoprotein_mean);
  load_double("oncoprotein_standard_deviation", oncoprotein_standard_deviation);
  load_double("oxygen_saturation_for_proliferation",
              oxygen_saturation_for_proliferation);
  load_double("oxygen_limit_for_proliferation", oxygen_limit_for_proliferation);
  load_double("oxygen_limit_for_necrosis", oxygen_limit_for_necrosis);
  load_double("oxygen_limit_for_necrosis_maximum",
              oxygen_limit_for_necrosis_maximum);
  load_double("time_lysis", time_lysis);
  load_double("maximum_necrosis_rate", maximum_necrosis_rate);

  load_double("default_oxygen_consumption_tumor_cell",
              default_oxygen_consumption_tumor_cell);
  load_double("default_volume_new_tumor_cell", default_volume_new_tumor_cell);
  load_double("default_volume_nucleus_tumor_cell",
              default_volume_nucleus_tumor_cell);
  load_double("default_fraction_fluid_tumor_cell",
              default_fraction_fluid_tumor_cell);
  load_double("average_time_transformation_random_rate",
              average_time_transformation_random_rate);
  load_double("standard_deviation_transformation_random_rate",
              standard_deviation_transformation_random_rate);
  load_double("adhesion_time", adhesion_time);
  load_double("oncoprotein_limit", oncoprotein_limit);
  load_double("oncoprotein_saturation", oncoprotein_saturation);

  // Difference between saturation and limit. This is always calculated here
  oncoprotein_difference = oncoprotein_saturation - oncoprotein_limit;

  load_double("volume_relaxation_rate_alive_tumor_cell_cytoplasm",
              volume_relaxation_rate_alive_tumor_cell_cytoplasm);
  load_double("volume_relaxation_rate_alive_tumor_cell_nucleus",
              volume_relaxation_rate_alive_tumor_cell_nucleus);
  load_double("volume_relaxation_rate_alive_tumor_cell_fluid",
              volume_relaxation_rate_alive_tumor_cell_fluid);

  load_double("volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell",
              volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell);
  load_double("volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell",
              volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell);
  load_double("volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell",
              volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell);

  load_double("volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell",
              volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell);
  load_double("volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell",
              volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell);
  load_double("volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell",
              volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell);

  load_double("threshold_cancer_cell_type1", threshold_cancer_cell_type1);
  load_double("threshold_cancer_cell_type2", threshold_cancer_cell_type2);
  load_double("threshold_cancer_cell_type3", threshold_cancer_cell_type3);
  load_double("threshold_cancer_cell_type4", threshold_cancer_cell_type4);

  if (jfile.contains("average_maximum_time_untill_apoptosis_cart")) {
    average_maximum_time_untill_apoptosis_cart =
        jfile["average_maximum_time_untill_apoptosis_cart"].get<double>();
  } else {
    average_maximum_time_untill_apoptosis_cart =
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
        dt_cycle * 10.0 * 24.0 * 60.0 / (dt_cycle + 1.0);
  }

  load_double("default_oxygen_consumption_cart",
              default_oxygen_consumption_cart);
  load_double("default_volume_new_cart_cell", default_volume_new_cart_cell);
  load_double("kill_rate_cart", kill_rate_cart);
  load_double("adhesion_rate_cart", adhesion_rate_cart);
  load_double("max_adhesion_distance_cart", max_adhesion_distance_cart);
  load_double("min_adhesion_distance_cart", min_adhesion_distance_cart);
  load_double("minimum_distance_from_tumor_to_spawn_cart",
              minimum_distance_from_tumor_to_spawn_cart);
  load_double("persistence_time_cart", persistence_time_cart);
  load_double("migration_bias_cart", migration_bias_cart);
  load_double("migration_speed_cart", migration_speed_cart);
  load_double("elastic_constant_cart", elastic_constant_cart);

  //
  // Computed constants that should not be directly changed
  //
  // Calculate steps per cycle. This is always calculated here
  steps_per_cell_cycle = static_cast<int>(dt_cycle / dt_step);
  // Calculate steps per day. This is always calculated here
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  steps_in_one_day = static_cast<size_t>(24 * 60 / dt_step);
  // Calculate the volume of a single mechanical voxel in μm³
  voxel_volume =
      (static_cast<real_t>(bounded_space_length) / resolution_grid_substances) *
      (static_cast<real_t>(bounded_space_length) / resolution_grid_substances) *
      (static_cast<real_t>(bounded_space_length) / resolution_grid_substances);
  // 1-migration_bias_cart
  migration_one_minus_bias_cart = 1.0 - migration_bias_cart;
  // Probability of a CAR-T cell to migrate in a given
  // mechanical time step
  motility_probability_cart = dt_mechanics / persistence_time_cart;
  // Probability of a Tumor cell to escape in a given
  // mechanical time step
  probability_escape_from_cart =
      dt_mechanics / (adhesion_time + kEpsilonProbability);
  // Maximum adhesion distance squared
  squared_max_adhesion_distance_cart =
      max_adhesion_distance_cart * max_adhesion_distance_cart;
  // Difference between min and max adhesion distance
  difference_cart_adhesion_distances =
      max_adhesion_distance_cart - min_adhesion_distance_cart;
  // Radius tumor cell
  radius_tumor_cell =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      std::cbrt(default_volume_new_tumor_cell * 3. / (4. * Math::kPi));
  // Radius cart cell
  radius_cart_cell =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      std::cbrt(default_volume_new_cart_cell * 3. / (4. * Math::kPi));
  // Max Distance for considering two cells as neighbours for force calculations
  // in μm (twice cell radius times max_relative_adhesion_distance + 0.1 to
  // avoid mismatch because of numerical errors)**2
  squared_max_distance_neighbors_force =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      std::pow(0.1 + 2 * radius_tumor_cell * max_relative_adhesion_distance, 2);

  //
  // Last constant that is derived from others but can be changed directly
  //
  // maximum squared distance to avoid CAR-T pushing
  // tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
  // CAR-T cell will only move to the tumor cell with the adhesion forces
  // (radiusCART + radiusTumorCell + 1 to avoid numerical errors)**2
  max_squared_distance_cart_moving_towards_tumor_cell =
      std::pow(radius_cart_cell + radius_tumor_cell + 1, 2);
}

void SimParam::PrintParams() const {
  std::cout << "\n========================================\n";
  std::cout << "        SIMULATION PARAMETERS\n";
  std::cout << "========================================\n\n";

  ///
  /// General Hyperparameters
  ///
  std::cout << "///\n";
  std::cout << "/// General Hyperparameters\n";
  std::cout << "///\n\n";

  std::cout << "Seed for random number generation: " << seed << "\n";
  std::cout << "Output Performance Statistics: "
            << (output_performance_statistics ? "true" : "false") << "\n";
  std::cout << "Total simulation time in minutes (30 days): "
            << total_minutes_to_simulate << "\n";
  std::cout << "Initial radius of the spherical tumor (micrometers): "
            << initial_tumor_radius << "\n";
  std::cout << "Length of the bounded space (micrometers): "
            << bounded_space_length << "\n\n";

  /// Treatment Dosages
  std::cout << "/// Treatment Dosages\n";
  std::cout << "///\n\n";
  std::cout << "CAR-T cell infusion schedule:\n";
  if (treatment.empty()) {
    std::cout << "  No treatment scheduled\n";
  } else {
    for (const auto& [day, dose] : treatment) {
      std::cout << "  Day " << day << ": " << dose << " CAR-T cells\n";
    }
  }
  std::cout << "\n";

  /// Time steps
  std::cout << "/// Time steps\n";
  std::cout << "///\n\n";
  std::cout << "Time step for substances secretion/consumption (minutes): "
            << dt_substances << "\n";
  std::cout << "Time step for the cell mechanics (minutes): " << dt_mechanics
            << "\n";
  std::cout << "Time step for the cell cycle (minutes): " << dt_cycle << "\n";
  std::cout << "General time step for the simulation: " << dt_step << "\n";
  std::cout << "Output CSV interval: " << output_csv_interval << "\n\n";

  /// Apoptotic cells volume change
  std::cout << "/// Apoptotic cells volume change\n";
  std::cout << "///\n\n";
  std::cout << "Volume relaxation rate cytoplasm apoptotic: "
            << volume_relaxation_rate_cytoplasm_apoptotic_cells << "\n";
  std::cout << "Volume relaxation rate nucleus apoptotic: "
            << volume_relaxation_rate_nucleus_apoptotic_cells << "\n";
  std::cout << "Volume relaxation rate fluid apoptotic: "
            << volume_relaxation_rate_fluid_apoptotic_cells << "\n";
  std::cout << "Time in minutes until an apoptotic cell is removed: "
            << time_apoptosis << "\n";
  std::cout << "Reduction of consumption rate of dead cells when they enter "
               "necrosis: "
            << reduction_consumption_dead_cells << "\n\n";

  /// Chemicals
  ///
  std::cout << "/// Chemicals\n";
  std::cout << "///\n\n";
  std::cout << "Number of voxels per axis for the substances grid: "
            << resolution_grid_substances << "\n";
  std::cout << "Diffusion coefficient of oxygen (μm²/min): "
            << diffusion_coefficient_oxygen << "\n";
  std::cout << "Decay constant of oxygen (min⁻¹): " << decay_constant_oxygen
            << "\n";
  std::cout << "Diffusion coefficient of immunostimulatory factor (μm²/min): "
            << diffusion_coefficient_immunostimulatory_factor << "\n";
  std::cout << "Decay constant of immunostimulatory factor (min⁻¹): "
            << decay_constant_immunostimulatory_factor << "\n";
  std::cout << "Reference level of oxygen at the boundaries (mmHg): "
            << oxygen_reference_level << "\n";
  std::cout << "Initial oxygen concentration in each voxel (mmHg): "
            << initial_oxygen_level << "\n";
  std::cout << "Oxygen saturation in the microenvironment (mmHg): "
            << oxygen_saturation << "\n\n";

  /// Forces
  ///
  std::cout << "/// Forces\n";
  std::cout << "///\n\n";
  std::cout << "Repulsion coefficient between tumor cells: "
            << cell_repulsion_between_tumor_tumor << "\n";
  std::cout << "Repulsion coefficient between CAR-T cells: "
            << cell_repulsion_between_cart_cart << "\n";
  std::cout << "Repulsion coefficient between CAR-T cells and tumor cells: "
            << cell_repulsion_between_cart_tumor << "\n";
  std::cout << "Repulsion coefficient between tumor cells and CAR-T cells: "
            << cell_repulsion_between_tumor_cart << "\n";
  std::cout << "Maximum relative adhesion distance for cell interactions: "
            << max_relative_adhesion_distance << "\n";
  std::cout << "Adhesion coefficient between tumor cells: "
            << cell_adhesion_between_tumor_tumor << "\n";
  std::cout << "Adhesion coefficient between CAR-T cells: "
            << cell_adhesion_between_cart_cart << "\n";
  std::cout << "Adhesion coefficient between CAR-T cells and tumor cells: "
            << cell_adhesion_between_cart_tumor << "\n";
  std::cout << "Adhesion coefficient between tumor cells and CAR-T cells: "
            << cell_adhesion_between_tumor_cart << "\n";
  std::cout << "Box length for mechanics calculations (micrometers): "
            << length_box_mechanics << "\n";
  std::cout
      << "Coefficient for the current time step in Adams-Bashforth method: "
      << dnew << "\n";
  std::cout
      << "Coefficient for the previous time step in Adams-Bashforth method: "
      << dold << "\n";
  std::cout << "Maximum squared distance to avoid CAR-T pushing tumor cells: "
            << max_squared_distance_cart_moving_towards_tumor_cell << "\n\n";

  ///
  /// TumorCell Hyperparameters
  ///
  std::cout << "///\n";
  std::cout << "/// TumorCell Hyperparameters\n";
  std::cout << "///\n\n";

  std::cout << "Rate of secretion of immunostimulatory factor of tumor cells "
               "per minute: "
            << rate_secretion_immunostimulatory_factor << "\n";
  std::cout
      << "Saturation density of immunostimulatory factor for tumor cells: "
      << saturation_density_immunostimulatory_factor << "\n";
  std::cout << "Mean level of oncoprotein expression in tumor cells: "
            << oncoprotein_mean << "\n";
  std::cout << "Standard deviation of oncoprotein expression in tumor cells: "
            << oncoprotein_standard_deviation << "\n";
  std::cout << "Oxygen saturation level in tumor cells for proliferation: "
            << oxygen_saturation_for_proliferation << "\n";
  std::cout << "Limit of oxygen level for tumor cell proliferation: "
            << oxygen_limit_for_proliferation << "\n";
  std::cout << "Limit of oxygen to start causing necrosis: "
            << oxygen_limit_for_necrosis << "\n";
  std::cout << "Limit of oxygen to maximum necrosis probability: "
            << oxygen_limit_for_necrosis_maximum << "\n";
  std::cout << "Time in minutes until a lysed necrotic cell is removed: "
            << time_lysis << "\n";
  std::cout << "Maximum rate per minute of necrosis for tumor cells: "
            << maximum_necrosis_rate << "\n";
  std::cout << "Default oxygen consumption rate of tumor cell: "
            << default_oxygen_consumption_tumor_cell << "\n";

  std::cout << "\nVolume parameters:\n";
  std::cout << "Default total volume of a new tumor cell (μm³): "
            << default_volume_new_tumor_cell << "\n";
  std::cout << "Default volume of the nucleus of a new tumor cell (μm³): "
            << default_volume_nucleus_tumor_cell << "\n";
  std::cout << "Default fraction of fluid volume in a new tumor cell: "
            << default_fraction_fluid_tumor_cell << "\n";

  std::cout << "\nTransformation and adhesion parameters:\n";
  std::cout << "Average time for transformation Random Rate (hours): "
            << average_time_transformation_random_rate << "\n";
  std::cout << "Standard Deviation for transformation Random Rate (hours): "
            << standard_deviation_transformation_random_rate << "\n";
  std::cout
      << "Average adhesion time for Tumor Cell under CAR-T attack (minutes): "
      << adhesion_time << "\n";
  std::cout << "Min oncoprotein level to be killed by a CAR-T cell: "
            << oncoprotein_limit << "\n";
  std::cout << "Max oncoprotein level: " << oncoprotein_saturation << "\n";
  std::cout << "Oncoprotein difference (saturation - limit): "
            << oncoprotein_difference << "\n";

  std::cout << "\nVolume relaxation rates for alive cells:\n";
  std::cout << "Volume relaxation rate alive cytoplasm (min⁻¹): "
            << volume_relaxation_rate_alive_tumor_cell_cytoplasm << "\n";
  std::cout << "Volume relaxation rate alive nucleus (min⁻¹): "
            << volume_relaxation_rate_alive_tumor_cell_nucleus << "\n";
  std::cout << "Volume relaxation rate alive fluid (min⁻¹): "
            << volume_relaxation_rate_alive_tumor_cell_fluid << "\n";

  std::cout << "\nVolume relaxation rates for necrotic swelling:\n";
  std::cout << "Volume relaxation rate cytoplasm necrotic swelling (min⁻¹): "
            << volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell
            << "\n";
  std::cout << "Volume relaxation rate nucleus necrotic swelling (min⁻¹): "
            << volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell
            << "\n";
  std::cout << "Volume relaxation rate fluid necrotic swelling (min⁻¹): "
            << volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell
            << "\n";

  std::cout << "\nVolume relaxation rates for necrotic lysed:\n";
  std::cout << "Volume relaxation rate cytoplasm necrotic lysed (min⁻¹): "
            << volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell
            << "\n";
  std::cout << "Volume relaxation rate nucleus necrotic lysed (min⁻¹): "
            << volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell << "\n";
  std::cout << "Volume relaxation rate fluid necrotic lysed (min⁻¹): "
            << volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell << "\n";

  std::cout << "\nThresholds in oncoprotein levels for differentiating cancer "
               "cell types:\n";
  std::cout << "Threshold Cancer Cell Type 1: " << threshold_cancer_cell_type1
            << "\n";
  std::cout << "Threshold Cancer Cell Type 2: " << threshold_cancer_cell_type2
            << "\n";
  std::cout << "Threshold Cancer Cell Type 3: " << threshold_cancer_cell_type3
            << "\n";
  std::cout << "Threshold Cancer Cell Type 4: " << threshold_cancer_cell_type4
            << "\n\n";

  ///
  /// CAR-T Cell Hyperparameters
  ///
  std::cout << "///\n";
  std::cout << "/// CAR-T Cell Hyperparameters\n";
  std::cout << "///\n\n";

  std::cout << "Average time in minutes until a CAR-T cell dies: "
            << average_maximum_time_untill_apoptosis_cart << "\n";
  std::cout << "Default oxygen consumption rate of CAR-T cell: "
            << default_oxygen_consumption_cart << "\n";

  std::cout << "\nVolume parameters:\n";
  std::cout << "Default total volume of a new CAR-T cell (μm³): "
            << default_volume_new_cart_cell << "\n";

  std::cout << "\nKilling and adhesion rates:\n";
  std::cout << "How often a CAR-T cell tries to kill an attached cancer cell "
               "(1/min): "
            << kill_rate_cart << "\n";
  std::cout
      << "How often a CAR-T cell tries to attach to a cancer cell (1/min): "
      << adhesion_rate_cart << "\n";
  std::cout << "Maximum adhesion distance between CAR-T and tumor cells "
               "(micrometers): "
            << max_adhesion_distance_cart << "\n";
  std::cout << "Minimum adhesion distance between CAR-T and tumor cells "
               "(micrometers): "
            << min_adhesion_distance_cart << "\n";

  std::cout << "\nSpawning parameters:\n";
  std::cout << "Minimum distance from the tumor for spawning CAR-T cells "
               "(micrometers): "
            << minimum_distance_from_tumor_to_spawn_cart << "\n";

  std::cout << "\nMotility parameters:\n";
  std::cout << "Average persistence time before CAR-T cell moves: "
            << persistence_time_cart << "\n";
  std::cout << "Migration bias (higher values = more directed movement): "
            << migration_bias_cart << "\n";
  std::cout << "Migration speed: " << migration_speed_cart << "\n";
  std::cout << "Elastic constant: " << elastic_constant_cart << "\n\n";

  std::cout << "========================================\n";
  std::cout << "      END OF PARAMETERS SUMMARY\n";
  std::cout << "========================================\n\n";
}

}  // namespace bdm
