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

#include "utils_aux.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "cart_cell.h"
#include "core/agent/agent.h"
#include "core/container/math_array.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <ios>
#include <tuple>
#include <vector>

namespace bdm {

// Samples a Gaussian value with given mean and standard deviation but all
// negative values are mapped to zero
real_t SamplePositiveGaussian(float mean, float sigma) {
  auto* random = Simulation::GetActive()->GetRandom();
  real_t value = random->Gaus(mean, sigma);
  if (value < 0.) {
    value = 0.;
  }
  return value;
}

std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius) {
  // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
  const real_t cell_radius =
      std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi) / kHalf;

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

//Function to compute the number of tumor cells of each type, the radius of the tumor, the average oncoprotein level and the average oxygen level of the cancer cells.
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, size_t, real_t, real_t, real_t> AnalyzeTumor() {

  auto* rm = Simulation::GetActive()->GetResourceManager();
  auto* oxygen_dgrid = rm->GetDiffusionGrid("oxygen");

  size_t total_num_tumor_cells = 0;
  size_t num_tumor_cells_type1 = 0;
  size_t num_tumor_cells_type2 = 0;
  size_t num_tumor_cells_type3 = 0;
  size_t num_tumor_cells_type4 = 0;
  size_t num_tumor_cells_type5_dead = 0;
  size_t num_alive_cart = 0;

  real_t max_dist_sq = 0.0;
  real_t acumulator_oncoprotein=0.0;
  real_t acumulator_oxygen_cancer_cells=0.0;

  rm->ForEachAgent([&](const Agent* agent) {
    if (auto* tumor_cell = dynamic_cast<const TumorCell*>(agent)) {
      total_num_tumor_cells++;
      const auto& pos = agent->GetPosition();

      // Accumulate oxygen level for average calculation
      acumulator_oxygen_cancer_cells += oxygen_dgrid->GetValue(pos);

      //computing tumor radius
      real_t dist_sq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
      if (dist_sq > max_dist_sq) {
        max_dist_sq = dist_sq;
      }

      // Count tumor cells by type and accumulate oncoprotein levels if the cell is alive
      switch (tumor_cell->GetType()) {
        case TumorCellType::kType1: acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel(); num_tumor_cells_type1++; break;
        case TumorCellType::kType2: acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel(); num_tumor_cells_type2++; break;
        case TumorCellType::kType3: acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel(); num_tumor_cells_type3++; break;
        case TumorCellType::kType4: acumulator_oncoprotein += tumor_cell->GetOncoproteinLevel(); num_tumor_cells_type4++; break;
        case TumorCellType::kType5: num_tumor_cells_type5_dead++; break;
        default: break;
      }
    } else if (auto* cart_cell = dynamic_cast<const CarTCell*>(agent)) {
      if (cart_cell->GetState() == CarTCellState::kAlive)
        num_alive_cart++;
    }
  });
  int total_alive_cells = num_tumor_cells_type1 + num_tumor_cells_type2 + num_tumor_cells_type3 + num_tumor_cells_type4;
  real_t average_oncoprotein = (total_alive_cells > 0) ? (acumulator_oncoprotein / total_alive_cells) : 0.0;
  real_t average_oxygen_cancer_cells = (total_num_tumor_cells > 0) ? (acumulator_oxygen_cancer_cells / total_num_tumor_cells) : 0.0;
  return {total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, num_alive_cart,std::sqrt(max_dist_sq), average_oncoprotein, average_oxygen_cancer_cells};
}

// Function to output summary CSV
void OutputSummary::operator()() {

  auto* simulation = Simulation::GetActive();
  auto* scheduler = simulation->GetScheduler();
  auto current_step = scheduler->GetSimulatedSteps();

  if (current_step % frequency_ == 0) {

    std::ofstream file("output/final_data.csv", std::ios::app);
    if (file.is_open()) {
      if (current_step == 0) {
        file << "total_days,total_hours,total_minutes,tumor_radius,num_cells,num_tumor_cells,tumor_cells_type1,tumor_cells_type2,tumor_cells_type3,tumor_cells_type4,tumor_cells_type5_dead,num_alive_cart,average_oncoprotein,average_oxygen_cancer_cells\n";// Header for CSV file
      }

      // Calculate time in days, hours, minutes
      double total_minutes = current_step*kDt;
      double total_hours = total_minutes / 60.0;
      double total_days = total_hours / 24.0;

      // Count total cells, tumor cells of each type and tumor radius
      size_t total_num_tumor_cells;
      size_t num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead,num_alive_cart;
      real_t tumor_radius, average_oncoprotein, average_oxygen_cancer_cells;
      std::tie(total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, num_alive_cart, tumor_radius, average_oncoprotein, average_oxygen_cancer_cells) = AnalyzeTumor();
      real_t total_num_cells=simulation->GetResourceManager()->GetNumAgents();

      //If a dosage is administred this exact time the numbers are not seen in 
      //the resource manager yet because of how BioDynaMo is built
      //therefore we need to add the just added new cells to the statistics here.
      size_t current_day_int = static_cast<size_t>(total_days);
      if (current_step % kStepsOneDay==0 && kTreatment.find(current_day_int) != kTreatment.end()){
        size_t just_spawned_cells= kTreatment.at(current_day_int);
        total_num_cells += just_spawned_cells;
        num_alive_cart += just_spawned_cells;
      }
      

      // Write data to CSV file
      file << total_days << ","
       << total_hours << ","
       << total_minutes << ","
       << tumor_radius << ","
       << total_num_cells << ","
       << total_num_tumor_cells << ","
       << num_tumor_cells_type1 << ","
       << num_tumor_cells_type2 << ","
       << num_tumor_cells_type3 << ","
       << num_tumor_cells_type4 << ","
       << num_tumor_cells_type5_dead << ","
       << num_alive_cart << ","
       << average_oncoprotein << ","
       << average_oxygen_cancer_cells << "\n";
    }
  }
}

// Function to spawn CAR-T cell dosages
void SpawnCart::operator()() {
  auto* simulation = Simulation::GetActive();
  auto* scheduler = simulation->GetScheduler();
  auto current_step = scheduler->GetSimulatedSteps();
  //This function only executes each day
  if (current_step % frequency_ != 0)
    return;
  //See if there is any dosage to apply in this day
  size_t current_day = current_step*kDt/ (60.0*24);
  size_t cells_to_spawn=0;

  if (kTreatment.find(current_day) != kTreatment.end())
    cells_to_spawn = kTreatment.at(current_day);

  //if there are cells to spawn in the treatment
  if (cells_to_spawn > 0) {

    //compute tumor radius
    auto* rm = simulation->GetResourceManager();
    real_t max_dist_sq = 0.0;

    rm->ForEachAgent([&](const Agent* agent) {
      if (std::strcmp(agent->GetTypeName(), "TumorCell") == 0) {
        const auto& pos = agent->GetPosition();
        real_t dist_sq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
        if (dist_sq > max_dist_sq) {
          max_dist_sq = dist_sq;
        }
      }
    });

    //the car-t spawns at least kMinimumDistanceCarTFromTumor micrometers away from the tumor
    real_t minimum_squared_radius= std::sqrt(max_dist_sq)+kMinimumDistanceCarTFromTumor;
    minimum_squared_radius *=minimum_squared_radius;

    //for generating car-t positions
    auto* rng= simulation->GetRandom();
    auto* param = simulation->GetParam();
    real_t min_b= param->min_bound;
    real_t max_b= param->max_bound;
    real_t px,py,pz,radi_sq;

    auto* ctxt = simulation->GetExecutionContext();
    CarTCell* cart = nullptr;

    for( unsigned int i=0 ;i < cells_to_spawn ; i++ ){
      //look for a valid position
      do{
          px = rng->Uniform(min_b, max_b);
          py = rng->Uniform(min_b, max_b);
          pz = rng->Uniform(min_b, max_b);
          radi_sq = px*px + py*py + pz*pz;
      }while(minimum_squared_radius > radi_sq);
      //spawn CAR-T
      cart = new CarTCell({px,py,pz});
      cart->AddBehavior(new StateControlCart());
      ctxt->AddAgent(cart);
    }
	}

}

//Function to generate a random direction unitary vector
Real3 GenerateRandomDirection() {
  auto* rnd= Simulation::GetActive()->GetRandom();
  // Generate a random direction vector
  double theta = kTwicePi * rnd->Uniform(0.0, 1.0);
  double phi = Math::kPi * rnd->Uniform(0.0, 1.0);

  double sin_phi = std::sin(phi);
  double cos_phi = std::cos(phi);

  return {sin_phi*std::cos(theta), sin_phi*std::sin(theta), cos_phi};
}

}  // namespace bdm
