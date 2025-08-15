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
      std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi) / 2;

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

// Function to compute the number of tumor cells of each type and the radius of
// the tumor
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, real_t>
ComputeNumberTumorCellsAndRadius() {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  size_t total_num_tumor_cells = 0;
  size_t num_tumor_cells_type1 = 0;
  size_t num_tumor_cells_type2 = 0;
  size_t num_tumor_cells_type3 = 0;
  size_t num_tumor_cells_type4 = 0;
  size_t num_tumor_cells_type5_dead = 0;

  real_t max_dist_sq = 0.0;

  rm->ForEachAgent([&](const Agent* agent) {
    if (const auto* tumor_cell = dynamic_cast<const TumorCell*>(agent)) {
      total_num_tumor_cells++;
      const auto& pos = agent->GetPosition();
      const real_t dist_sq =
          pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
      if (dist_sq > max_dist_sq) {
        max_dist_sq = dist_sq;
      }

      // Count tumor cells by type
      switch (tumor_cell->GetType()) {
        case TumorCellType::kType1:
          num_tumor_cells_type1++;
          break;
        case TumorCellType::kType2:
          num_tumor_cells_type2++;
          break;
        case TumorCellType::kType3:
          num_tumor_cells_type3++;
          break;
        case TumorCellType::kType4:
          num_tumor_cells_type4++;
          break;
        case TumorCellType::kType5:
          num_tumor_cells_type5_dead++;
          break;
        default:
          break;
      }
    }
  });
  return {total_num_tumor_cells, num_tumor_cells_type1,
          num_tumor_cells_type2, num_tumor_cells_type3,
          num_tumor_cells_type4, num_tumor_cells_type5_dead,
          std::sqrt(max_dist_sq)};
}

// Function to output summary CSV
void OutputSummary::operator()() {
  auto* scheduler = Simulation::GetActive()->GetScheduler();
  auto current_step = scheduler->GetSimulatedSteps();

  if (current_step % frequency_ == 0) {
    std::ofstream file("output/final_data.csv", std::ios::app);
    if (file.is_open()) {
      // Header for the CSV file
      if (current_step == 0) {
        file << "total_days,total_hours,total_minutes,tumor_radius,num_cells,"
                "num_tumor_cells,tumor_cells_type1,tumor_cells_type2,tumor_"
                "cells_type3,tumor_cells_type4,tumor_cells_type5_dead\n";
      }

      // Calculate time in days, hours, minutes
      const double total_minutes =
          Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
      const double total_hours = total_minutes / kMinutesInAnHour;
      const double total_days = total_hours / kHoursInADay;

      // Count total cells, tumor cells of each type and tumor radius
      size_t total_num_tumor_cells = 0;
      size_t num_tumor_cells_type1 = 0;
      size_t num_tumor_cells_type2 = 0;
      size_t num_tumor_cells_type3 = 0;
      size_t num_tumor_cells_type4 = 0;
      size_t num_tumor_cells_type5_dead = 0;
      real_t tumor_radius = 0.0;

      std::tie(total_num_tumor_cells, num_tumor_cells_type1,
               num_tumor_cells_type2, num_tumor_cells_type3,
               num_tumor_cells_type4, num_tumor_cells_type5_dead,
               tumor_radius) = ComputeNumberTumorCellsAndRadius();

      // Write data to CSV file
      file << total_days << "," << total_hours << "," << total_minutes << ","
           << tumor_radius << ","
           << Simulation::GetActive()->GetResourceManager()->GetNumAgents()
           << "," << total_num_tumor_cells << "," << num_tumor_cells_type1
           << "," << num_tumor_cells_type2 << "," << num_tumor_cells_type3
           << "," << num_tumor_cells_type4 << "," << num_tumor_cells_type5_dead
           << "\n";
    }
  }
}

}  // namespace bdm
