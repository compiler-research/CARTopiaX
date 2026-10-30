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

#include "params/hyperparams.h"
#include "core/param/param.h"
#include "core/param/param_group.h"
#include "core/real_t.h"
#include "core/util/log.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

namespace bdm {

// Define the static member kUid for SimParam
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

// Function to read a JSON file and load the parameters
// The parameters that are not found in the file keep their default
void SimParam::LoadParams(const std::string& filename) {
  std::ifstream file(filename);
  std::stringstream content;
  // If the json can be opened, read it, otherwise keep all the defaults
  if (file.is_open()) {
    content << file.rdbuf();
  }
  LoadParamsFromString(content.str());
}

// Function to load the parameters from a JSON string
// The parameters that are not found in the string keep their default
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void SimParam::LoadParamsFromString(const std::string& json_string) {
  nlohmann::json jfile;
  // If the string is not empty, try to parse it, otherwise leave jfile empty
  if (!json_string.empty()) {
    try {
      jfile = nlohmann::json::parse(json_string);
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

  auto load_string = [&](const std::string& key, std::string& var) {
    if (jfile.contains(key)) {
      var = jfile[key].get<std::string>();
    }
  };

  // Load parameters from JSON file
  load_int("seed", seed);
  load_bool("output_performance_statistics", output_performance_statistics);
  load_int("total_minutes_to_simulate", total_minutes_to_simulate);

  std::string tumor_shape_string;
  load_string("tumor_shape", tumor_shape_string);
  if (jfile.contains("tumor_shape")) {
    tumor_shape_string = jfile["tumor_shape"].get<std::string>();
    if (tumor_shape_string == "sphere") {
      tumor_shape = TumorShape::kSphere;
    } else if (tumor_shape_string == "cylinder") {
      tumor_shape = TumorShape::kCylinder;
    } else {
      tumor_shape = TumorShape::kSphere;
      Log::Error("SimParam::LoadParams",
                 "Unknown tumor shape '" + tumor_shape_string +
                     "'. It should be either 'sphere' or 'cylinder'. Loading "
                     "default value ('sphere') instead");
    }
  }

  load_int("initial_number_of_cylindrical_tumor_cells",
           initial_number_of_cylindrical_tumor_cells);
  load_int("bounded_space_length", bounded_space_length);
  std::string bound_space_toplogy_string;
  if (jfile.contains("bound_space_toplogy")) {
    bound_space_toplogy_string =
        jfile["bound_space_toplogy"].get<std::string>();
    if (bound_space_toplogy_string == "torus") {
      bound_space_toplogy = Param::BoundSpaceMode::kTorus;
    } else if (bound_space_toplogy_string == "open") {
      bound_space_toplogy = Param::BoundSpaceMode::kOpen;
    } else if (bound_space_toplogy_string == "closed") {
      bound_space_toplogy = Param::BoundSpaceMode::kClosed;
    } else {
      bound_space_toplogy = Param::BoundSpaceMode::kTorus;
      Log::Error("SimParam::LoadParams",
                 "Unknown bound space topology '" + bound_space_toplogy_string +
                     "'. It should be either 'torus', 'open' or 'closed'. "
                     "Loading default value ('torus') instead");
    }
  }

  load_double("initial_spherical_tumor_radius", initial_spherical_tumor_radius);
  load_double("cylindrical_tumor_radius", cylindrical_tumor_radius);
  load_double("cylindrical_tumor_height", cylindrical_tumor_height);

  // Half of the bounded space length (integer division, as in the original
  // defaults)
  const int half_bounded_space_length = bounded_space_length / 2;

  if (jfile.contains("bounded_space_min_allowed_z")) {
    bounded_space_min_allowed_z =
        jfile["bounded_space_min_allowed_z"].get<double>();
  } else {
    bounded_space_min_allowed_z = -half_bounded_space_length;
  }
  if (jfile.contains("bounded_space_max_allowed_z")) {
    bounded_space_max_allowed_z =
        jfile["bounded_space_max_allowed_z"].get<double>();
  } else {
    bounded_space_max_allowed_z = half_bounded_space_length;
  }

  if (jfile.contains("bounded_space_max_allowed_radius")) {
    bounded_space_max_allowed_radius =
        jfile["bounded_space_max_allowed_radius"].get<double>();
  } else {
    // if it is not defined we assign the whole length of the bounded space to
    // avoid any restrictions
    bounded_space_max_allowed_radius = bounded_space_length;
  }

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

  load_bool("output_information_dependent_on_radius",
            output_information_dependent_on_radius);

  load_double("max_radius_analysis_csv_dependent_on_radius",
              max_radius_analysis_csv_dependent_on_radius);

  load_int("num_radius_intervals", num_radius_intervals);

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

  if (jfile.contains("min_initial_z_substances")) {
    min_initial_z_substances = jfile["min_initial_z_substances"].get<double>();
  } else {
    min_initial_z_substances = -half_bounded_space_length;
  }

  if (jfile.contains("max_initial_z_substances")) {
    max_initial_z_substances = jfile["max_initial_z_substances"].get<double>();
  } else {
    max_initial_z_substances = half_bounded_space_length;
  }

  if (jfile.contains("lateral_oxygen_production_min_z")) {
    lateral_oxygen_production_min_z =
        jfile["lateral_oxygen_production_min_z"].get<double>();
  } else {
    lateral_oxygen_production_min_z = -half_bounded_space_length;
  }

  if (jfile.contains("lateral_oxygen_production_max_z")) {
    lateral_oxygen_production_max_z =
        jfile["lateral_oxygen_production_max_z"].get<double>();
  } else {
    lateral_oxygen_production_max_z = half_bounded_space_length;
  }
  load_bool("diffuse_oxygen_on_z_axis", diffuse_oxygen_on_z_axis);
  load_bool("diffuse_immunostimulatory_factor_on_z_axis",
            diffuse_immunostimulatory_factor_on_z_axis);
  load_bool("diffuse_glucose_on_z_axis", diffuse_glucose_on_z_axis);

  load_double("diffusion_coefficient_oxygen", diffusion_coefficient_oxygen);
  load_double("decay_constant_oxygen", decay_constant_oxygen);
  load_double("oxygen_reference_level", oxygen_reference_level);
  load_double("initial_oxygen_level", initial_oxygen_level);
  load_bool("add_immunostimulatory_factor", add_immunostimulatory_factor);
  load_double("diffusion_coefficient_immunostimulatory_factor",
              diffusion_coefficient_immunostimulatory_factor);
  load_double("decay_constant_immunostimulatory_factor",
              decay_constant_immunostimulatory_factor);
  load_bool("add_glucose", add_glucose);
  load_double("diffusion_coefficient_glucose", diffusion_coefficient_glucose);
  load_double("decay_constant_glucose", decay_constant_glucose);
  load_double("initial_glucose_level", initial_glucose_level);

  if (jfile.contains("max_radius_glucose_initialization")) {
    max_radius_glucose_initialization =
        jfile["max_radius_glucose_initialization"].get<double>();
  } else {
    max_radius_glucose_initialization = bounded_space_length;
  }

  if (jfile.contains("diffuse_glucose_on_z_axis")) {
    diffuse_glucose_on_z_axis = jfile["diffuse_glucose_on_z_axis"].get<bool>();
  } else {
    // if the tumor shape is cylindrical it should be set to false, otherwise it
    // should be set to true
    diffuse_glucose_on_z_axis = (tumor_shape != TumorShape::kCylinder);
  }

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
  load_double("glucose_saturation_for_tumor_cell_growth",
              glucose_saturation_for_tumor_cell_growth);
  load_double("glucose_limit_for_tumor_cell_growth",
              glucose_limit_for_tumor_cell_growth);
  load_double("oxygen_limit_for_necrosis", oxygen_limit_for_necrosis);
  load_double("oxygen_limit_for_necrosis_maximum",
              oxygen_limit_for_necrosis_maximum);
  load_double("maximum_necrosis_lack_of_oxygen_rate",
              maximum_necrosis_lack_of_oxygen_rate);
  load_double("glucose_limit_for_death", glucose_limit_for_death);
  load_double("glucose_limit_for_death_maximum",
              glucose_limit_for_death_maximum);
  load_double("maximum_death_lack_of_glucose_rate",
              maximum_death_lack_of_glucose_rate);
  load_double("basal_death_probability_cancer_cells",
              basal_death_probability_cancer_cells);
  load_double("time_lysis", time_lysis);

  load_double("default_oxygen_consumption_tumor_cell",
              default_oxygen_consumption_tumor_cell);
  load_double("default_glucose_consumption_tumor_cell",
              default_glucose_consumption_tumor_cell);
  load_double("default_volume_new_tumor_cell", default_volume_new_tumor_cell);
  load_double("std_volume_new_tumor_cell", std_volume_new_tumor_cell);
  load_double("max_volume_new_tumor_cell", max_volume_new_tumor_cell);
  load_double("min_volume_new_tumor_cell", min_volume_new_tumor_cell);
  load_double("default_fraction_of_volume_for_nucleus_tumor_cell",
              default_fraction_of_volume_for_nucleus_tumor_cell);
  load_double("default_fraction_fluid_tumor_cell",
              default_fraction_fluid_tumor_cell);
  load_double("average_time_transformation_random_rate",
              average_time_transformation_random_rate);
  load_double("standard_deviation_transformation_random_rate",
              standard_deviation_transformation_random_rate);
  load_double("minimum_tumor_cell_target_volume_fraction_for_division",
              minimum_tumor_cell_target_volume_fraction_for_division);

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

  if (jfile.contains("average_maximum_time_until_apoptosis_cart")) {
    average_maximum_time_until_apoptosis_cart =
        jfile["average_maximum_time_until_apoptosis_cart"].get<double>();
  } else {
    average_maximum_time_until_apoptosis_cart =
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
        dt_cycle * 10.0 * 24.0 * 60.0 / (dt_cycle + 1.0);
  }

  load_double("default_oxygen_consumption_cart",
              default_oxygen_consumption_cart);
  load_double("default_glucose_consumption_cart",
              default_glucose_consumption_cart);
  load_double("default_volume_new_cart_cell", default_volume_new_cart_cell);
  load_double("std_volume_new_cart_cell", std_volume_new_cart_cell);
  load_double("max_volume_new_cart_cell", max_volume_new_cart_cell);
  load_double("min_volume_new_cart_cell", min_volume_new_cart_cell);
  load_double("default_fraction_fluid_cart_cell",
              default_fraction_fluid_cart_cell);
  load_double("default_fraction_of_volume_for_nucleus_cart_cell",
              default_fraction_of_volume_for_nucleus_cart_cell);
  load_double("kill_rate_cart", kill_rate_cart);
  load_double("adhesion_rate_cart", adhesion_rate_cart);
  load_double("max_adhesion_distance_cart", max_adhesion_distance_cart);
  load_double("min_adhesion_distance_cart", min_adhesion_distance_cart);
  load_double("minimum_distance_from_tumor_to_spawn_cart",
              minimum_distance_from_tumor_to_spawn_cart);
  load_double("persistence_time_cart", persistence_time_cart);
  load_double("avg_migration_bias_cart", avg_migration_bias_cart);
  load_double("std_migration_bias_cart", std_migration_bias_cart);
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
  // Calculate the volume of a single substance voxel in μm³
  voxel_volume =
      (static_cast<real_t>(bounded_space_length) / resolution_grid_substances) *
      (static_cast<real_t>(bounded_space_length) / resolution_grid_substances) *
      (static_cast<real_t>(bounded_space_length) / resolution_grid_substances);

  // Squared max allowed radius
  bounded_space_max_allowed_radius_squared =
      bounded_space_max_allowed_radius * bounded_space_max_allowed_radius;

  // Probability of a CAR-T cell to migrate in a given
  // mechanical time step
  motility_probability_cart =
      (persistence_time_cart == 0) ? 1.0 : dt_mechanics / persistence_time_cart;
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
  // maximum squared distance to limit CAR-T direct attraction to tumor cells
  // tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
  // CAR-T cell will only move to the tumor cell with the adhesion forces
  // By default (radiusCART + radiusTumorCell + 1 to avoid numerical errors)**2,
  // unless it is given in the JSON file
  if (jfile.contains("max_squared_distance_cart_moving_towards_tumor_cell")) {
    max_squared_distance_cart_moving_towards_tumor_cell =
        jfile["max_squared_distance_cart_moving_towards_tumor_cell"]
            .get<double>();
  } else {
    max_squared_distance_cart_moving_towards_tumor_cell =
        std::pow(radius_cart_cell + radius_tumor_cell + 1, 2);
  }
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
  switch (tumor_shape) {
    case TumorShape::kSphere:
      std::cout << "Tumor shape: 'sphere' \n";
      std::cout << "Initial radius of the spherical tumor (micrometers): "
                << initial_spherical_tumor_radius << "\n";
      break;
    case TumorShape::kCylinder:
      std::cout << "Tumor shape: 'cylinder' \n";
      std::cout << "Initial radius of the cylindrical tumor (micrometers): "
                << cylindrical_tumor_radius << "\n";
      std::cout << "Initial height of the cylindrical tumor (micrometers): "
                << cylindrical_tumor_height << "\n";
      std::cout << "Min allowed Z coordinate for cells (micrometers): "
                << bounded_space_min_allowed_z << "\n";
      std::cout << "Max allowed Z coordinate for cells (micrometers): "
                << bounded_space_max_allowed_z << "\n";
      std::cout << "Max allowed radius for cells (micrometers): "
                << bounded_space_max_allowed_radius << "\n";
      std::cout << "Initial number of cylindrical tumor cells: "
                << initial_number_of_cylindrical_tumor_cells << "\n";
      break;
  }
  std::cout << "Length of the bounded space (micrometers): "
            << bounded_space_length << "\n\n";
  switch (bound_space_toplogy) {
    case Param::BoundSpaceMode::kTorus:
      std::cout << "Bound space topology: 'torus' \n";
      break;
    case Param::BoundSpaceMode::kOpen:
      std::cout << "Bounded space topology: 'open' \n";
      break;
    case Param::BoundSpaceMode::kClosed:
      std::cout << "Bounded space topology: 'closed' \n";
      break;
  }

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
  std::cout << "Output information dependent on radius: "
            << (output_information_dependent_on_radius ? "true" : "false")
            << "\n\n";
  if (output_information_dependent_on_radius) {
    std::cout << "Number of radius intervals to output the agregated "
                 "information CSV: "
              << num_radius_intervals << "\n\n";
    std::cout << "Maximum radius for the analysis of information dependent on "
                 "radius: "
              << max_radius_analysis_csv_dependent_on_radius << "\n\n";
  }

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
  std::cout << "Reduction of consumption rate of dead cells when they die: "
            << reduction_consumption_dead_cells << "\n\n";

  /// Chemicals
  ///
  std::cout << "/// Chemicals\n";
  std::cout << "///\n\n";
  std::cout << "Number of voxels per axis for the substances grid: "
            << resolution_grid_substances << "\n";
  std::cout << "Minimum initial z-coordinate for substances values different "
               "than 0 (micrometers): "
            << min_initial_z_substances << "\n";
  std::cout << "Maximum initial z-coordinate for substances values different "
               "than 0 (micrometers): "
            << max_initial_z_substances << "\n";
  std::cout
      << "Minimum z-coordinate for lateral oxygen production (micrometers): "
      << lateral_oxygen_production_min_z << "\n";
  std::cout
      << "Maximum z-coordinate for lateral oxygen production (micrometers): "
      << lateral_oxygen_production_max_z << "\n";
  std::cout << "Whether to diffuse oxygen on the z-axis: "
            << (diffuse_oxygen_on_z_axis ? "true" : "false") << "\n";
  std::cout << "Whether to diffuse immunostimulatory factor on the z-axis: "
            << (diffuse_immunostimulatory_factor_on_z_axis ? "true" : "false")
            << "\n";
  std::cout << "Whether to diffuse glucose on the z-axis: "
            << (diffuse_glucose_on_z_axis ? "true" : "false") << "\n";
  std::cout << "Diffusion coefficient of oxygen (μm²/min): "
            << diffusion_coefficient_oxygen << "\n";
  std::cout << "Decay constant of oxygen (min⁻¹): " << decay_constant_oxygen
            << "\n";
  std::cout << "Reference level of oxygen at the boundaries (mmHg): "
            << oxygen_reference_level << "\n";
  std::cout << "Initial oxygen concentration in each voxel (mmHg): "
            << initial_oxygen_level << "\n";
  std::cout << "Whether to add immunostimulatory factor: "
            << (add_immunostimulatory_factor ? "true" : "false") << "\n";
  std::cout << "Diffusion coefficient of immunostimulatory factor (μm²/min): "
            << diffusion_coefficient_immunostimulatory_factor << "\n";
  std::cout << "Decay constant of immunostimulatory factor (min⁻¹): "
            << decay_constant_immunostimulatory_factor << "\n";
  std::cout << "Whether to add glucose: " << (add_glucose ? "true" : "false")
            << "\n";
  std::cout << "Diffusion coefficient of glucose (μm²/min): "
            << diffusion_coefficient_glucose << "\n";
  std::cout << "Decay constant of glucose (min⁻¹): " << decay_constant_glucose
            << "\n";
  std::cout << "Initial glucose concentration in each voxel (mM): "
            << initial_glucose_level << "\n\n";
  std::cout << "Maximum radius for glucose initialization (micrometers): "
            << max_radius_glucose_initialization << "\n\n";

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
  std::cout << "Glucose saturation level in tumor cells for growth: "
            << glucose_saturation_for_tumor_cell_growth << "\n";
  std::cout << "Limit of glucose level for tumor cell growth: "
            << glucose_limit_for_tumor_cell_growth << "\n";
  std::cout << "Limit of oxygen to start causing necrosis: "
            << oxygen_limit_for_necrosis << "\n";
  std::cout << "Limit of oxygen to maximum necrosis probability: "
            << oxygen_limit_for_necrosis_maximum << "\n";
  std::cout << "Basal death probability for tumor cells: "
            << basal_death_probability_cancer_cells << "\n";
  std::cout << "Time in minutes until a lysed necrotic cell is removed: "
            << time_lysis << "\n";
  std::cout << "Maximum rate per minute of necrosis for tumor cells: "
            << maximum_necrosis_lack_of_oxygen_rate << "\n";
  std::cout << "Limit of glucose to start causing death: "
            << glucose_limit_for_death << "\n";
  std::cout << "Limit of glucose to maximum death probability: "
            << glucose_limit_for_death_maximum << "\n";
  std::cout << "Maximum rate per minute of death for tumor cells: "
            << maximum_death_lack_of_glucose_rate << "\n";

  std::cout << "Default oxygen consumption rate of tumor cell: "
            << default_oxygen_consumption_tumor_cell << "\n";
  std::cout << "Default glucose consumption rate of tumor cell: "
            << default_glucose_consumption_tumor_cell << "\n";

  std::cout << "\nVolume parameters:\n";
  std::cout << "Default total volume of a new tumor cell (μm³): "
            << default_volume_new_tumor_cell << "\n";
  std::cout
      << "Standard deviation of the total volume of a new tumor cell (μm³): "
      << std_volume_new_tumor_cell << "\n";
  std::cout << "Maximum total volume of a new tumor cell (μm³): "
            << max_volume_new_tumor_cell << "\n";
  std::cout << "Minimum total volume of a new tumor cell (μm³): "
            << min_volume_new_tumor_cell << "\n";
  std::cout << "Default fraction of volume for nucleus tumor cell: "
            << default_fraction_of_volume_for_nucleus_tumor_cell << "\n";
  std::cout << "Default fraction of fluid volume in a new tumor cell: "
            << default_fraction_fluid_tumor_cell << "\n";

  std::cout << "\nTransformation and adhesion parameters:\n";
  std::cout << "Average time for transformation Random Rate (hours): "
            << average_time_transformation_random_rate << "\n";
  std::cout << "Standard Deviation for transformation Random Rate (hours): "
            << standard_deviation_transformation_random_rate << "\n";
  std::cout << "Minimum tumor cell target volume fraction for division: "
            << minimum_tumor_cell_target_volume_fraction_for_division << "\n";
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
            << average_maximum_time_until_apoptosis_cart << "\n";
  std::cout << "Default oxygen consumption rate of CAR-T cell: "
            << default_oxygen_consumption_cart << "\n";
  std::cout << "Default glucose consumption rate of CAR-T cell: "
            << default_glucose_consumption_cart << "\n";

  std::cout << "\nVolume parameters:\n";
  std::cout << "Default total volume of a new CAR-T cell (μm³): "
            << default_volume_new_cart_cell << "\n";
  std::cout
      << "Standard deviation of the total volume of a new CAR-T cell (μm³): "
      << std_volume_new_cart_cell << "\n";
  std::cout << "Maximum total volume of a new CAR-T cell (μm³): "
            << max_volume_new_cart_cell << "\n";
  std::cout << "Minimum total volume of a new CAR-T cell (μm³): "
            << min_volume_new_cart_cell << "\n";
  std::cout << "Default fraction of fluid volume in a new CAR-T cell: "
            << default_fraction_fluid_cart_cell << "\n";
  std::cout << "Default fraction of volume for nucleus CAR-T cell: "
            << default_fraction_of_volume_for_nucleus_cart_cell << "\n";

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
  std::cout
      << "Average Migration bias (higher values = more directed movement): "
      << avg_migration_bias_cart << "\n";
  std::cout << "Standard deviation migration bias: " << std_migration_bias_cart
            << "\n";
  std::cout << "Migration speed: " << migration_speed_cart << "\n";
  std::cout << "Elastic constant: " << elastic_constant_cart << "\n\n";

  std::cout << "========================================\n";
  std::cout << "      END OF PARAMETERS SUMMARY\n";
  std::cout << "========================================\n\n";
}

}  // namespace bdm
