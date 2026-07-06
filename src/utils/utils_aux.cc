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

#include "utils/utils_aux.h"
#include "agents/cart_cell.h"
#include "agents/tumor_cell.h"
#include "params/hyperparams.h"
#include "core/agent/agent.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/param/param.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include "core/scheduler.h"
#include "core/simulation.h"
#include "core/util/math.h"
#include "core/util/random.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <memory>
#include <tuple>
#include <vector>

namespace bdm {

// Samples a Gaussian value with given mean and standard deviation but all
// negative values are mapped to zero
real_t SamplePositiveGaussian(real_t mean, real_t sigma) {
  Random* random = Simulation::GetActive()->GetRandom();
  real_t value = random->Gaus(mean, sigma);
  if (value < 0.) {
    value = 0.;
  }
  return value;
}

std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius) {
  const auto* sparams = Simulation::GetActive()->GetParam()->Get<SimParam>();
  // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
  const real_t cell_radius =
      std::cbrt(sparams->default_volume_new_tumor_cell * 6 / Math::kPi) / kHalf;

  std::vector<Real3> positions;

  // Hexagonal close-packing spacing
  const real_t spacing_x = cell_radius * std::sqrt(3.0);
  const real_t spacing_y = cell_radius * 2.0;
  const real_t spacing_z = cell_radius * std::sqrt(3.0);

  // Use integer counters instead of floating-point loop variables
  const int z_steps = static_cast<int>((2 * sphere_radius) / spacing_z) + 1;
  const int x_steps = static_cast<int>((2 * sphere_radius) / spacing_x) + 1;
  const int y_steps = static_cast<int>((2 * sphere_radius) / spacing_y) + 1;

  for (int zi = 0; zi < z_steps; ++zi) {
    const real_t z = -sphere_radius + zi * spacing_z;
    for (int xi = 0; xi < x_steps; ++xi) {
      const real_t x = -sphere_radius + xi * spacing_x;
      for (int yi = 0; yi < y_steps; ++yi) {
        const real_t y = -sphere_radius + yi * spacing_y;

        // Compute cell center with HCP offset
        const real_t px = x + (zi % 2) * 0.5 * cell_radius;
        const real_t py = y + (xi % 2) * cell_radius;
        const real_t pz = z;

        const real_t dist = std::sqrt(px * px + py * py + pz * pz);

        if (dist <= sphere_radius) {
          positions.push_back({px, py, pz});
        }
      }
    }
  }

  return positions;
}

/// Create a cylindrical arrangement of tumor cells randomly distributed within its volume
std::vector<Real3> CreateCylinderOfTumorCells(real_t cylinder_radius, real_t cylinder_height, size_t number_of_cells) {
  
  std::vector<Real3> positions;
  positions.reserve(number_of_cells);

  Random* random = Simulation::GetActive()->GetRandom();

  for (size_t i = 0; i < number_of_cells; ++i) {
    // sqrt ensures uniform distribution over the disk area (area element = r dr dθ)
    real_t angle  = random->Uniform(0.0, kTwicePi);
    real_t radius = cylinder_radius * std::sqrt(random->Uniform(0.0, 1.0));
    real_t height = random->Uniform(-cylinder_height / 2.0, cylinder_height / 2.0);

    positions.push_back({radius * std::cos(angle),
                         radius * std::sin(angle),
                         height});
  }

  return positions;
}

// Function to compute the number of tumor cells of each type, the radius of the
// tumor, the number of alive and dead cart cells, the average oncoprotein level, the average oxygen level of the
// cancer cells and the average oxygen level of all cells. Only considering cells within the specified inner and outer radius.
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, real_t,
           real_t, real_t, real_t>
AnalyzeTumor(real_t inner_radius_considered, real_t outer_radius_considered) {
  real_t inner_radius_squared = inner_radius_considered * inner_radius_considered;
  real_t outer_radius_squared = outer_radius_considered * outer_radius_considered;
  
  Simulation* sim = Simulation::GetActive();
  ResourceManager* rm = sim->GetResourceManager();
  DiffusionGrid* oxygen_dgrid = rm->GetDiffusionGrid("oxygen");
  const auto* sparams = sim->GetParam()->Get<SimParam>(); 

  int total_num_tumor_cells = 0;
  int num_tumor_cells_type1 = 0;
  int num_tumor_cells_type2 = 0;
  int num_tumor_cells_type3 = 0;
  int num_tumor_cells_type4 = 0;
  int num_tumor_cells_type5_dead = 0;
  int num_alive_cart = 0;
  int num_dead_cart = 0;

  real_t max_dist_sq = 0.0;
  real_t acumulator_oncoprotein = 0.0;
  real_t acumulator_oxygen_cancer_cells = 0.0;
  real_t acumulator_oxygen_all_cells = 0.0;

  bool tumor_shape_is_cylindrical = sparams->tumor_shape == "cylinder";
  bool tumor_shape_is_spherical = sparams->tumor_shape == "sphere";

  rm->ForEachAgent([&](const Agent* agent) {
    // Compute the distance to the center depending on the tumor shape
    const Real3& pos = agent->GetPosition();
    real_t dist_sq = 0.0;
    if (tumor_shape_is_cylindrical) {
      // Only consider x and y for cylindrical rumor, distance to the axis of the cylinder
      dist_sq = pos[0] * pos[0] + pos[1] * pos[1];  
    } else if (tumor_shape_is_spherical) {
      // Consider all three dimensions for spherical radius
      dist_sq = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];  
    } else{
      Log::Error("AnalyzeTumor", "Unknown tumor shape, please use 'sphere' or 'cylinder'.");
    } 

    //If the cell is outside the considered radius range, skip it
    if (dist_sq < inner_radius_squared || dist_sq >= outer_radius_squared) {
      return;
    }

    // The agent is within the considered radius range, analyze it
    if (const auto* tumor_cell = dynamic_cast<const TumorCell*>(agent)) {

      // Accumulate oxygen level for average calculation
      acumulator_oxygen_all_cells += oxygen_dgrid->GetValue(pos);

      total_num_tumor_cells++;

      // Accumulate oxygen level for average calculation
      acumulator_oxygen_cancer_cells += oxygen_dgrid->GetValue(pos);


      // computing tumor radius
      if (dist_sq > max_dist_sq) {
        max_dist_sq = dist_sq;
      }

      // Count tumor cells by type and accumulate oncoprotein levels if the cell
      // is alive
      switch (tumor_cell->GetType()) {
        case TumorCellType::kType1:
          acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel();
          num_tumor_cells_type1++;
          break;
        case TumorCellType::kType2:
          acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel();
          num_tumor_cells_type2++;
          break;
        case TumorCellType::kType3:
          acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel();
          num_tumor_cells_type3++;
          break;
        case TumorCellType::kType4:
          acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel();
          num_tumor_cells_type4++;
          break;
        case TumorCellType::kType5:
          num_tumor_cells_type5_dead++;
          break;
        default:
          break;
      }
    } else if (const auto* cart_cell = dynamic_cast<const CarTCell*>(agent)) {
      if (cart_cell->GetState() == CarTCellState::kAlive) {
        // Accumulate oxygen level for average calculation
        acumulator_oxygen_all_cells += oxygen_dgrid->GetValue(pos);
        num_alive_cart++;
      } else if (cart_cell->GetState() == CarTCellState::kApoptotic) {
        // Accumulate oxygen level for average calculation
        acumulator_oxygen_all_cells += oxygen_dgrid->GetValue(pos);
        //Dead CART cell
        num_dead_cart++;
      }
    }
  });
  const int total_alive_cells = num_tumor_cells_type1 + num_tumor_cells_type2 +
                                num_tumor_cells_type3 + num_tumor_cells_type4;
  const real_t average_oncoprotein =
      (total_alive_cells > 0) ? (acumulator_oncoprotein / total_alive_cells)
                              : 0.0;
  const real_t average_oxygen_cancer_cells =
      (total_num_tumor_cells > 0)
          ? (acumulator_oxygen_cancer_cells / total_num_tumor_cells)
          : 0.0;
  const real_t average_oxygen_all_cells =
      (total_num_tumor_cells + num_alive_cart + num_dead_cart > 0)
          ? (acumulator_oxygen_all_cells / (total_num_tumor_cells + num_alive_cart + num_dead_cart))
          : 0.0;
  return {total_num_tumor_cells, num_tumor_cells_type1,
          num_tumor_cells_type2, num_tumor_cells_type3,
          num_tumor_cells_type4, num_tumor_cells_type5_dead,
          num_alive_cart, num_dead_cart, std::sqrt(max_dist_sq),
          average_oncoprotein,   average_oxygen_cancer_cells, average_oxygen_all_cells};
}

// Function to output summary CSV
void OutputSummary::operator()() {
  Simulation* simulation = Simulation::GetActive();
  const auto* sparams = simulation->GetParam()->Get<SimParam>();
  Scheduler* scheduler = simulation->GetScheduler();
  const uint64_t current_step = scheduler->GetSimulatedSteps();

  if (current_step % frequency_ == 0) {
    // Calculate time in days, hours, minutes
    const double total_minutes = static_cast<double>(current_step) * sparams->dt_step;
    const double total_hours = total_minutes / kMinutesInAnHour;
    const double total_days = total_hours / kHoursInADay;

    // Delete csv content current_step == 0 to, otherwise append mode
    std::ofstream file("output/final_data.csv",
                       current_step == 0 ? std::ios::trunc : std::ios::app);
    if (file.is_open()) {
      if (current_step == 0) {
        file
            << "total_days,total_hours,total_minutes,tumor_radius,num_cells,"
               "num_tumor_cells,tumor_cells_type1,tumor_cells_type2,tumor_"
               "cells_type3,tumor_cells_type4,tumor_cells_type5_dead,num_alive_"
               "cart,num_dead_cart,average_oncoprotein,average_oxygen_cancer_cells,average_oxygen_all_cells\n";  // Header
                                                                          // for
                                                                          // CSV
                                                                          // file
      }

      // Count total cells, tumor cells of each type and tumor radius
      int total_num_tumor_cells = 0;
      int num_tumor_cells_type1 = 0;
      int num_tumor_cells_type2 = 0;
      int num_tumor_cells_type3 = 0;
      int num_tumor_cells_type4 = 0;
      int num_tumor_cells_type5_dead = 0;
      int num_alive_cart = 0;
      int num_dead_cart = 0;
      real_t tumor_radius = 0.0;
      real_t average_oncoprotein = 0.0;
      real_t average_oxygen_cancer_cells = 0.0;
      real_t average_oxygen_all_cells = 0.0;
      // Analyze tumor with no limits on inner and outer radius
      std::tie(total_num_tumor_cells, num_tumor_cells_type1,
               num_tumor_cells_type2, num_tumor_cells_type3,
               num_tumor_cells_type4, num_tumor_cells_type5_dead,
               num_alive_cart, num_dead_cart, tumor_radius, average_oncoprotein,
               average_oxygen_cancer_cells, average_oxygen_all_cells) = AnalyzeTumor(0, sparams->bounded_space_length);
      size_t total_num_cells = simulation->GetResourceManager()->GetNumAgents();

      // If a dosage is administred this exact time the numbers are not seen in
      // the resource manager yet because of how BioDynaMo is built
      // therefore we need to add the just added new cells to the statistics
      // here.
      const auto current_day_int = static_cast<int>(total_days);
      if (current_step % sparams->steps_in_one_day == 0 &&
          sparams->treatment.find(current_day_int) !=
              sparams->treatment.end()) {
        const size_t just_spawned_cells =
            sparams->treatment.at(current_day_int);
        total_num_cells += just_spawned_cells;
        num_alive_cart += static_cast<int>(just_spawned_cells);
      }

      // Write data to CSV file
      file << total_days << "," << total_hours << "," << total_minutes << ","
           << tumor_radius << "," << total_num_cells << ","
           << total_num_tumor_cells << "," << num_tumor_cells_type1 << ","
           << num_tumor_cells_type2 << "," << num_tumor_cells_type3 << ","
           << num_tumor_cells_type4 << "," << num_tumor_cells_type5_dead << ","
           << num_alive_cart << ","<< num_dead_cart << "," << average_oncoprotein << ","
           << average_oxygen_cancer_cells << "," << average_oxygen_all_cells << "\n";
    }

    if (sparams->output_information_dependent_on_radius) {
      OutputInformationBasedOnRadiusCSV(current_step, total_minutes, total_hours, total_days);
    }
  }
}

// Function to output information based on radius to CSV
void OutputInformationBasedOnRadiusCSV(const uint64_t current_step, const real_t total_minutes, const real_t total_hours, const real_t total_days) {
  const Simulation* simulation = Simulation::GetActive();
  const auto* sparams = simulation->GetParam()->Get<SimParam>();
  const real_t interval_size = sparams->max_radius_analysis_csv_dependent_on_radius / sparams->num_radius_intervals;

  // Delete csv content current_step == 0 to, otherwise append mode
  std::ofstream file("output/data_dependent_on_radius_tumor.csv",
                      current_step == 0 ? std::ios::trunc : std::ios::app);
  if (file.is_open()) {
    if (current_step == 0) {
      // Header for CSV file
      file << "total_days,total_hours,total_minutes";
      for (int i = 0; i < sparams->num_radius_intervals; ++i) {
        file << ",num_alive_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",num_dead_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",total_num_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",num_alive_tumor_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",total_num_tumor_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",num_alive_cart_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",num_dead_cart_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",tumor_cells_type1_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",tumor_cells_type2_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",tumor_cells_type3_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",tumor_cells_type4_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",tumor_cells_type5_dead_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",average_oncoprotein_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",average_oxygen_cancer_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;
        file << ",average_oxygen_all_cells_radius_" << i * interval_size << "_to_" << (i + 1) * interval_size;  
      }
      // End of header line
      file << "\n";
    }
    // Write the time information for the current step
    file << total_days << "," << total_hours << "," << total_minutes;

    for (int i = 0; i < sparams->num_radius_intervals; ++i) {


      real_t inner_radius_considered = i * interval_size;
      real_t outer_radius_considered = (i + 1) * interval_size;

      // Count total cells, tumor cells of each type and tumor radius
      int total_num_tumor_cells = 0;
      int num_tumor_cells_type1 = 0;
      int num_tumor_cells_type2 = 0;
      int num_tumor_cells_type3 = 0;
      int num_tumor_cells_type4 = 0;
      int num_tumor_cells_type5_dead = 0;
      int num_alive_cart = 0;
      int num_dead_cart = 0;
      real_t tumor_radius = 0.0;
      real_t average_oncoprotein = 0.0;
      real_t average_oxygen_cancer_cells = 0.0;
      real_t average_oxygen_all_cells = 0.0;
      // Analyze tumor with no limits on inner and outer radius
      std::tie(total_num_tumor_cells, num_tumor_cells_type1,
              num_tumor_cells_type2, num_tumor_cells_type3,
              num_tumor_cells_type4, num_tumor_cells_type5_dead,
              num_alive_cart, num_dead_cart, tumor_radius, average_oncoprotein,
              average_oxygen_cancer_cells, average_oxygen_all_cells) = AnalyzeTumor(inner_radius_considered, outer_radius_considered);

      const int num_alive_tumor_cells = num_tumor_cells_type1 + num_tumor_cells_type2 + num_tumor_cells_type3 + num_tumor_cells_type4;
      const int num_alive_cells = num_alive_cart + num_alive_tumor_cells;
      const int num_dead_cells = num_dead_cart + num_tumor_cells_type5_dead;
      const int total_num_cells = num_alive_cells + num_dead_cells;

      // Write data to CSV file
      file << "," << num_alive_cells
            << "," << num_dead_cells 
            << "," << total_num_cells
            << "," << num_alive_tumor_cells
            << "," << total_num_tumor_cells
            << "," << num_alive_cart
            << "," << num_dead_cart
            << "," << num_tumor_cells_type1
            << "," << num_tumor_cells_type2
            << "," << num_tumor_cells_type3
            << "," << num_tumor_cells_type4
            << "," << num_tumor_cells_type5_dead
            << "," << average_oncoprotein
            << "," << average_oxygen_cancer_cells
            << "," << average_oxygen_all_cells;
    }  
    // End of line for the current step
    file << "\n";
  }
}; 

// Function to spawn CAR-T cell dosages
void SpawnCart::operator()() {
  Simulation* simulation = Simulation::GetActive();
  const auto* sparams = simulation->GetParam()->Get<SimParam>();

  Scheduler* scheduler = simulation->GetScheduler();
  const uint64_t current_step = scheduler->GetSimulatedSteps();
  // This function only executes each day
  if (current_step % frequency_ != 0) {
    return;
  }
  // See if there is any dosage to apply in this day
  const auto current_day_int =
      static_cast<int>(static_cast<double>(current_step) * sparams->dt_step /
                       (kMinutesInAnHour * kHoursInADay));
  size_t cells_to_spawn = 0;

  if (sparams->treatment.find(current_day_int) != sparams->treatment.end()) {
    cells_to_spawn = sparams->treatment.at(current_day_int);
  }

  // if there are cells to spawn in the treatment
  if (cells_to_spawn > 0) {
    // compute tumor radius
    ResourceManager* rm = simulation->GetResourceManager();
    real_t max_dist_sq = 0.0;

    rm->ForEachAgent([&](const Agent* agent) {
      if (const auto* cancer_cell = dynamic_cast<const TumorCell*>(agent)) {
        const Real3& pos = cancer_cell->GetPosition();
        const real_t dist_sq =
            pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        if (dist_sq > max_dist_sq) {
          max_dist_sq = dist_sq;
        }
      }
    });

    // the car-t spawns at least
    // sparams->minimum_distance_from_tumor_to_spawn_cart micrometers away from
    // the tumor
    real_t minimum_squared_radius =
        std::sqrt(max_dist_sq) +
        sparams->minimum_distance_from_tumor_to_spawn_cart;
    minimum_squared_radius *= minimum_squared_radius;

    // for generating car-t positions
    Random* rng = simulation->GetRandom();
    const Param* param = simulation->GetParam();
    const real_t min_b = param->min_bound;
    const real_t max_b = param->max_bound;
    real_t px = 0.0;
    real_t py = 0.0;
    real_t pz = 0.0;
    real_t radi_sq = 0.0;

    ExecutionContext* ctxt = simulation->GetExecutionContext();

    for (unsigned int i = 0; i < cells_to_spawn; i++) {
      // look for a valid position
      while (true) {
        px = rng->Uniform(min_b, max_b);
        py = rng->Uniform(min_b, max_b);
        pz = rng->Uniform(min_b, max_b);
        radi_sq = px * px + py * py + pz * pz;
        if (radi_sq >= minimum_squared_radius) {
          break;
        }
      }
      // spawn CAR-T
      std::unique_ptr<CarTCell> cart =
          std::make_unique<CarTCell>(Real3{px, py, pz});
      std::unique_ptr<StateControlCart> state_control =
          std::make_unique<StateControlCart>();
      cart->AddBehavior(state_control.release());
      ctxt->AddAgent(cart.release());
    }
  }
}

// Function to generate a random direction unitary vector
Real3 GenerateRandomDirection() {
  Random* rnd = Simulation::GetActive()->GetRandom();
  // Generate a random direction vector
  const double theta = kTwicePi * rnd->Uniform(0.0, 1.0);
  const double phi = Math::kPi * rnd->Uniform(0.0, 1.0);

  const double sin_phi = std::sin(phi);
  const double cos_phi = std::cos(phi);

  return {sin_phi * std::cos(theta), sin_phi * std::sin(theta), cos_phi};
}

}  // namespace bdm
