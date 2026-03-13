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

#ifndef TUMOR_HYPERPARAMS_H_
#define TUMOR_HYPERPARAMS_H_

#include "core/param/param_group.h"
#include "core/real_t.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>
#include <map>
#include <string>

namespace bdm {

///
/// Some constants defined for code legibility
///
/// Twice Pi
constexpr real_t kTwicePi = 2. * Math::kPi;
/// Constant for division by 2.0
constexpr real_t kHalf = 2.0;
/// Epsilon for avoiding division by 0
constexpr real_t kEpsilon = 1e-10;
/// kEpsilon distance
constexpr real_t kEpsilonDistance = 1e-5;
/// Epsilon for avoiding division by 0 in probabilities
constexpr real_t kEpsilonProbability = 1e-15;
/// Large time to avoid division by 0
constexpr real_t kTimeTooLarge = 1e100;
/// Minutes in an hour
constexpr real_t kMinutesInAnHour = 60.0;
/// Hours in a day
constexpr real_t kHoursInADay = 24.0;

/// Contains the default values of the hyperparameters used in the simulation.
struct SimParam : public ParamGroup {
  // NOLINTNEXTLINE(modernize-type-traits,llvm-else-after-return,readability-else-after-return,cppcoreguidelines-owning-memory)
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  SimParam(const SimParam&) = default;
  SimParam(SimParam&&) = default;
  SimParam& operator=(const SimParam&) = default;
  SimParam& operator=(SimParam&&) = default;

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  ///
  /// General Hyperparameters
  ///

  /// Seed for random number generation
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int seed = 1;
  /// Output Performance Statistics
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  bool output_performance_statistics = false;
  /// Total simulation time in minutes (30 days by default)
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int total_minutes_to_simulate = 43200;
  /// Initial radius of the spherical tumor (group of cancer cells) in
  /// micrometers
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t initial_tumor_radius = 150;
  /// Length of the bounded space in micrometers
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int bounded_space_length = 1000;

  /// Treatment Dosages
  ///
  /// Specifies the CAR-T cell infusion schedule as a map where:
  ///   - The key represents the day of treatment (starting from day 0).
  ///   - The value represents the number of CAR-T cells administered on that
  ///   day.
  /// Example: On day 0 and day 8, 3957 CAR-T cells are introduced (matching the
  /// initial tumor cell count). {0, 3957}, {8, 3957}}
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  std::map<int, int> treatment = {{0, 3957}, {8, 3957}};

  /// Time steps
  /// Time step for substances secretion/consumption in minutes
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t dt_substances = 0.01;
  /// Time step for the cell mechanics in minutes
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t dt_mechanics = 0.1;
  /// Time step for the cell cycle in minutes
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t dt_cycle = 6;
  /// General time step for the simulation: it is the same as dt_mechanics
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t dt_step = 0.1;

  /// Output little summary each half output_csv_interval.
  /// By default half a day (12h * 60min / dt_step)
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int output_csv_interval = 7200;

  /// Apoptotic cells volume change
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_cytoplasm_apoptotic_cells = 0.0166667;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_nucleus_apoptotic_cells = 0.00583333;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_fluid_apoptotic_cells = 0.0;
  /// Time in minutes until an apoptotic cell is removed from the simulation
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t time_apoptosis = 516;
  /// Reduction of consumption rate of dead cells when they enter necrosis
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t reduction_consumption_dead_cells = 0.1;

  /// Chemicals
  ///  Number of voxels per axis for the substances grid
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int resolution_grid_substances = 50;
  /// Diffusion coefficient of oxygen in μm²/min
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t diffusion_coefficient_oxygen = 100000;
  /// Decay constant of oxygen in min⁻¹
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t decay_constant_oxygen = 0.1;
  /// Diffusion coefficient of immunostimulatory factor in μm²/min
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t diffusion_coefficient_immunostimulatory_factor = 1000;
  /// Decay constant of immunostimulatory factor in min⁻¹
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t decay_constant_immunostimulatory_factor = 0.016;
  /// Reference level of oxygen at the boundaries in mmHg
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oxygen_reference_level = 38;
  /// Initial oxygen concentration in each voxel in mmHg
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t initial_oxygen_level = 38;
  /// Oxygen saturation in the microenvironment in mmHg
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oxygen_saturation = 30;

  /// Forces
  ///  Repulsion coeficient between tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_repulsion_between_tumor_tumor = 10;
  /// Repulsion coeficient between CAR-T cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_repulsion_between_cart_cart = 50;
  /// Repulsion coeficient between CAR-T cells and tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_repulsion_between_cart_tumor = 50;
  /// Repulsion coeficient between tumor cells and CAR-T cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_repulsion_between_tumor_cart = 10;
  /// Maximum relative adhesion distance for cell interactions
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t max_relative_adhesion_distance = 1.25;
  /// Adhesion coeficient between tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_adhesion_between_tumor_tumor = 0.4;
  /// Adhesion coeficient between CAR-T cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_adhesion_between_cart_cart = 0;
  /// Adhesion coeficient between CAR-T cells and tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_adhesion_between_cart_tumor = 0;
  /// Adhesion coeficient between tumor cells and CAR-T cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t cell_adhesion_between_tumor_cart = 0;
  /// Box length for mechanics calculations (in micrometers). Smaller values
  /// improve simulation efficiency
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int length_box_mechanics = 22;
  // Coefficientes for the two step Adams-Bashforth approximation of the time
  // derivative for position position(t + dt) ≈ position(t) + dt * [ 1.5 *
  // velocity(t) - 0.5 * velocity(t - dt) ]
  /// Coefficient for the current time step in the Adams-Bashforth method (dt
  /// * 1.5)
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t dnew = 0.15;
  /// Coefficient for the previous time step in the Adams-Bashforth method (dt *
  /// -0.5)
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t dold = -0.05;
  /// maximum squared distance to avoid CAR-T pushing
  /// tumor cells If a CAR-T and a Tumor Cell are closer than this distance, the
  /// CAR-T cell will only move to the tumor cell with the adhesion forces
  /// (radiusCART + radiusTumorCell + 0.1 to avoid numerical errors)**2
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t max_squared_distance_cart_moving_towards_tumor_cell = 317.746;

  ///
  /// TumorCell Hyperparameters
  ///

  /// Rate of secretion of immunostimulatory factor of tumor cells per minute
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t rate_secretion_immunostimulatory_factor = 10;
  /// Saturation density of immunostimulatory factor for tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t saturation_density_immunostimulatory_factor = 1;
  /// Mean level of oncoprotein expression in tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oncoprotein_mean = 1;
  /// Standard deviation of oncoprotein expression in tumor cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oncoprotein_standard_deviation = 0.25;
  /// Oxygen saturation level in tumor cells for proliferation
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oxygen_saturation_for_proliferation = 38;
  /// Limit of oxygen level for tumor cell proliferation
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oxygen_limit_for_proliferation = 10;
  /// Limit of oxygen to start causing necrosis
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oxygen_limit_for_necrosis = 5;
  /// Limit of oxygen to maximum necrosis probability
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oxygen_limit_for_necrosis_maximum = 2.5;
  /// Time in minutes until a lysed necrotic cell is removed from the simulation
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t time_lysis = 86400;
  /// Maximum rate per minute of necrosis for tumor cells in case of hypoxia
  /// with 0 oxygen
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t maximum_necrosis_rate = 0.00277778;
  /// Default oxygen consumption rate of tumor cell
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t default_oxygen_consumption_tumor_cell = 10;
  /// Volume parameters
  /// Default total volume of a new tumor cell in μm³
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t default_volume_new_tumor_cell = 2494;
  /// Default volume of the nucleus of a new tumor cell in μm³
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t default_volume_nucleus_tumor_cell = 540;
  /// Default fraction of fluid volume in a new tumor cell
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t default_fraction_fluid_tumor_cell = 0.75;
  /// Average time for transformation Random Rate in hours
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t average_time_transformation_random_rate = 38.6;
  /// Standard Deviation for transformation Random Rate in hours
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t standard_deviation_transformation_random_rate = 3.7;
  /// Average adhesion time in minutes for Tumor Cell under CAR-T attack before
  /// escaping
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t adhesion_time = 60;
  /// Min oncoprotein level to be killed by a CAR-T cell
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oncoprotein_limit = 0.5;
  /// Max oncoprotein level
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oncoprotein_saturation = 2.0;
  /// Do not modify this line: difference between saturation and limit
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t oncoprotein_difference = 1.5;

  /// volume relaxation rate (min^-1) for each state
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_alive_tumor_cell_cytoplasm = 0.00216667;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_alive_tumor_cell_nucleus = 0.00366667;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_alive_tumor_cell_fluid = 0.0216667;

  real_t volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      5.33333e-05;
  real_t volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      0.000216667;
  real_t volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      0.000833333;

  real_t volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell =
      // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
      5.33333e-05;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell = 0.000216667;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell = 0.000833333;

  /// Thresholds in oncoprotein levels for differentiating 4 cancer cell types
  /// depending on their aggressiveness
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t threshold_cancer_cell_type1 = 1.5;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t threshold_cancer_cell_type2 = 1.0;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t threshold_cancer_cell_type3 = 0.5;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t threshold_cancer_cell_type4 = 0.0;

  ///
  /// CAR-T Cell Hyperparameters
  ///

  /// Average time in minutes until a CAR-T cell dies
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t average_maximum_time_untill_apoptosis_cart = 12342.86;
  /// Default oxygen consumption rate of CAR-T cell
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t default_oxygen_consumption_cart = 1;
  /// Volume parameters
  ///  Default total volume of a new CAR-T cell in μm³
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t default_volume_new_cart_cell = 2494;

  /// How often a CAR-T cell tries to kill an attached cancer cell in 1/min
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t kill_rate_cart = 0.06667;
  /// How often a CAR-T cell tries to attach to a cancer cell in 1/min
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t adhesion_rate_cart = 0.013;
  /// Maximum adhesion distance between CAR-T and tumor cells in micrometers
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t max_adhesion_distance_cart = 18;
  /// Minimum adhesion distance between CAR-T and tumor cells in micrometers
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t min_adhesion_distance_cart = 14;

  /// Minimum distance in micrometers from the tumor for spawning CAR-T cells
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t minimum_distance_from_tumor_to_spawn_cart = 50;

  /// Motility parameters
  /// Average persistence time before CAR-T cell moves
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t persistence_time_cart = 10;  // 10 minutes
  /// Higher bias (\in [0,1]) makes CAR-T movement more directed toward
  /// immunostimulatory factor source; while a bias of 0 makes the movement
  /// random
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t migration_bias_cart = 0.5;
  /// Migration speed
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t migration_speed_cart = 5;
  /// Elastic constant
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t elastic_constant_cart = 0.01;

  ///
  /// Fixed Constants that should not be directly changed
  ///

  /// Number of steps per cycle step Needs to be computed to avoid errors with
  /// fmod
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  int steps_per_cell_cycle = 60;
  /// Steps in one day
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  size_t steps_in_one_day = 14400;
  /// Volume of a single voxel in μm³
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t voxel_volume = 8000;
  /// 1-migration_bias_cart
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t migration_one_minus_bias_cart = 0.5;
  /// Probability of a CAR-T cell to migrate in a given
  /// mechanical time step
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t motility_probability_cart = 0.01;
  /// Probability of a Tumor cell to escape in a given
  /// mechanical time step
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t probability_escape_from_cart = 0.00166667;

  /// Maximum adhesion distance squared
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t squared_max_adhesion_distance_cart = 324;
  /// Difference between min and max adhesion distance
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t difference_cart_adhesion_distances = 4;

  /// Radius tumor cell
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t radius_tumor_cell = 8.41271;

  /// Radius cart cell
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t radius_cart_cell = 8.41271;

  /// Max Distance for considering two cells as neighbours for force
  /// calculations in μm (twice cell radius times max_relative_adhesion_distance
  /// + 0.1 to avoid mismatch because of numerical errors)**2
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  real_t squared_max_distance_neighbors_force = 446.552;
  // NOLINTEND(misc-non-private-member-variables-in-classes)

  ///
  /// Function to load the parameters
  ///

  /// Function to read a JSON file and load the hyperparameters in it.
  /// The parameters that are not found in the file will keep their
  /// default value
  void LoadParams(const std::string& filename);

  ///
  /// Function to print all the parameters with their values
  ///

  /// Print all the parameters with their values
  void PrintParams() const;
};

}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_
