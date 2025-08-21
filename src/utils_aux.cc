/*
 * Copyright 2025 compiler-research.org, Salvador de la Torre Gonzalez, Luciana Melina Luque
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

namespace bdm {

// Samples a Gaussian value with given mean and standard deviation but all negative values are mapped to zero
real_t SamplePositiveGaussian(float mean, float sigma) {
  auto* random = Simulation::GetActive()->GetRandom();
  real_t value = random->Gaus(mean, sigma);
  if(value < 0.) {value = 0.;}
  return value;
}


std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius) {
    // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
    real_t cell_radius = std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi)/2;

    std::vector<Real3> positions;

    // Hexagonal close-packing spacing
    real_t spacing_x = cell_radius * std::sqrt(3.0);
    real_t spacing_y = cell_radius * 2.0;
    real_t spacing_z = cell_radius * std::sqrt(3.0);

    int zc = 0;
    for (real_t z = -sphere_radius; z < sphere_radius; z += spacing_z, ++zc) {
        int xc = 0;
        for (real_t x = -sphere_radius; x < sphere_radius; x += spacing_x, ++xc) {
            int yc = 0;
            for (real_t y = -sphere_radius; y < sphere_radius; y += spacing_y, ++yc) {

                // Compute cell center with HCP offset
                real_t px = x + (zc % 2) * 0.5 * cell_radius;
                real_t py = y + (xc % 2) * cell_radius;
                real_t pz = z;

                real_t dist = std::sqrt(px * px + py * py + pz * pz);

                if (dist <= sphere_radius) {
                    positions.push_back({px, py, pz});
                }
            }
        }
    }

    return positions;
}

//Function to compute the number of tumor cells of each type and the radius of the tumor
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, real_t> ComputeNumberTumorCellsAndRadius() {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  size_t total_num_tumor_cells = 0;
  size_t num_tumor_cells_type1 = 0;
  size_t num_tumor_cells_type2 = 0;
  size_t num_tumor_cells_type3 = 0;
  size_t num_tumor_cells_type4 = 0;
  size_t num_tumor_cells_type5_dead = 0;

  real_t max_dist_sq = 0.0;

  rm->ForEachAgent([&](const Agent* agent) {
    if (auto* tumor_cell = dynamic_cast<const TumorCell*>(agent)) {
      total_num_tumor_cells++;
      const auto& pos = agent->GetPosition();
      real_t dist_sq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
      if (dist_sq > max_dist_sq) {
        max_dist_sq = dist_sq;
      }

      // Count tumor cells by type
      switch (tumor_cell->GetType()) {
        case 1: num_tumor_cells_type1++; break;
        case 2: num_tumor_cells_type2++; break;
        case 3: num_tumor_cells_type3++; break;
        case 4: num_tumor_cells_type4++; break;
        case 5: num_tumor_cells_type5_dead++; break;
        default: break;
      }
    }
  });
  return {total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, std::sqrt(max_dist_sq)};
}

// Function to output summary CSV
void OutputSummary::operator()() {
  auto* scheduler = Simulation::GetActive()->GetScheduler();
  auto current_step = scheduler->GetSimulatedSteps();

  if (current_step % frequency_ == 0) {
    std::ofstream file("output/final_data.csv", std::ios::app);
    if (file.is_open()) {
      if (current_step == 0) {
        file << "total_days,total_hours,total_minutes,tumor_radius,num_cells,num_tumor_cells,tumor_cells_type1,tumor_cells_type2,tumor_cells_type3,tumor_cells_type4,tumor_cells_type5_dead\n";// Header for CSV file
      }

      // Calculate time in days, hours, minutes
      double total_minutes = Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
      double total_hours = total_minutes / 60.0;
      double total_days = total_hours / 24.0;

      // Count total cells, tumor cells of each type and tumor radius
      size_t total_num_tumor_cells;
      size_t num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead;
      real_t tumor_radius;
      std::tie(total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, tumor_radius) = ComputeNumberTumorCellsAndRadius();
      // Write data to CSV file
      file << total_days << ","
       << total_hours << ","
       << total_minutes << ","
       << tumor_radius << ","
       << Simulation::GetActive()->GetResourceManager()->GetNumAgents() << ","//number of cells
       << total_num_tumor_cells << ","
       << num_tumor_cells_type1 << ","
       << num_tumor_cells_type2 << ","
       << num_tumor_cells_type3 << ","
       << num_tumor_cells_type4 << ","
       << num_tumor_cells_type5_dead << "\n";
    }
  }
}

}  // namespace bdm