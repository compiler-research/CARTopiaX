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

#include <algorithm>
#include <cmath>

#include "core/agent/agent.h"
#include "core/agent/cell.h"
#include "core/container/math_array.h"
#include "core/interaction_force.h"
#include "core/real_t.h"

#include "forces_tumor_cart.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "utils_aux.h"

namespace bdm {

Real4 InteractionVelocity::Calculate(const Agent* lhs, const Agent* rhs) const {
  const auto* a = dynamic_cast<const Cell*>(lhs);
  const auto* b = dynamic_cast<const Cell*>(rhs);

  // Ignore self-interaction
  if (a->GetUid() == b->GetUid()) {
    return {0.0, 0.0, 0.0, 0.0};
  }

  Real3 displacement = a->GetPosition() - b->GetPosition();

  // For periodic boundary conditions, we need to adjust the displacement
  displacement[0] =
      displacement[0] -
      (kBoundedSpaceLength)*round(displacement[0] / (kBoundedSpaceLength));
  displacement[1] =
      displacement[1] -
      (kBoundedSpaceLength)*round(displacement[1] / (kBoundedSpaceLength));
  displacement[2] =
      displacement[2] -
      (kBoundedSpaceLength)*round(displacement[2] / (kBoundedSpaceLength));

  const real_t dist_sq = displacement[0] * displacement[0] +
                         displacement[1] * displacement[1] +
                         displacement[2] * displacement[2];
  const real_t distance = std::max(std::sqrt(dist_sq), kEpsilonDistance);

  constexpr real_t kHalf = 2.0;
  const real_t radius_a = a->GetDiameter() / kHalf;
  const real_t radius_b = b->GetDiameter() / kHalf;
  const real_t combined_radius = radius_a + radius_b;
  // combined_radius=16.8254;//Debug
  // std::cout << "Debug: combined_radius = " << combined_radius << ", distance
  // = " << distance << std::endl;// Debug output
  real_t temp_r = 0.0;

  const auto* a_tumor = dynamic_cast<const TumorCell*>(a);
  const auto* b_tumor = dynamic_cast<const TumorCell*>(b);

  if (distance < combined_radius) {
    // 1 - d/combined_radius
    temp_r = 1.0 - distance / combined_radius;
    // (1 - d/combined_radius)^2
    temp_r *= temp_r;

    real_t repulsion = NAN;
    // std::cout << "temp_r = " << temp_r<< std::endl;// Debug output

    if ((a_tumor != nullptr) && (b_tumor != nullptr)) {  // two tumor cells
      repulsion = kRepulsionTumorTumor;  // std::sqrt(kRepulsionTumorTumor *
                                         // kRepulsionTumorTumor);
    } else if ((a_tumor == nullptr) &&
               (b_tumor == nullptr)) {  // two CAR-T cells
      repulsion =
          kRepulsionCartCart;  // std::sqrt(kRepulsionCartCart*kRepulsionCartCart);
    } else {                   // one tumor cell and one CAR-T
      repulsion = std::sqrt(kRepulsionCartTumor * kRepulsionTumorCart);
    }

    temp_r *= repulsion;
  }

  // std::cout << "temp_r after repulsion = " << temp_r<< std::endl;// Debug
  // output

  // Adhesion
  const real_t max_interaction_distance =
      kMaxRelativeAdhesionDistance * combined_radius;
  // max_interaction_distance=21.0318;//Debug
  // std::cout << "max_interaction_distance = " << max_interaction_distance <<
  // std::endl;// Debug output

  if (distance < max_interaction_distance) {
    // 1 - d/S
    real_t temp_a = 1.0 - distance / max_interaction_distance;
    // (1-d/S)^2
    temp_a *= temp_a;

    // std::cout << "temp_a = " << temp_a << std::endl;// Debug output

    real_t adhesion = NAN;                               // Initialize to NAN
    if ((a_tumor != nullptr) && (b_tumor != nullptr)) {  // two tumor cells
      adhesion = kAdhesionTumorTumor;
    } else if ((a_tumor == nullptr) &&
               (b_tumor == nullptr)) {  // two CAR-T cells
      adhesion = kAdhesionCartCart;
    } else {  // one tumor cell and one CAR-T
      adhesion = std::sqrt(kAdhesionCartTumor * kAdhesionTumorCart);
    }

    // std::cout << "adhesion = " << adhesion << std::endl;// Debug output

    temp_a *= adhesion;
    temp_r -= temp_a;

    // std::cout << "temp_a after adhesion= " << temp_a << std::endl;// Debug
    // output
  }

  if (std::abs(temp_r) < kEpsilon) {
    return {0.0, 0.0, 0.0, 0.0};
  }
  const real_t force_magnitude = temp_r / distance;

  // Debug Output volcities
  // std::ofstream file("output/intercation_velocities.csv", std::ios::app);
  // if (file.is_open()) {

  //   real_t total_minutes =
  //   Simulation::GetActive()->GetScheduler()->GetSimulatedTime(); Real3
  //   position = a->GetPosition();
  //   // Write data to CSV file
  //   file << total_minutes << ",position"
  //    << position[0] << ","
  //    << position[1] << ","
  //    << position[2] << ",displacement"
  //     << displacement[0] << ","
  //     << displacement[1] << ","
  //     << displacement[2] << ",distance"
  //     << distance << ",force_magnitude"
  //     << force_magnitude << ",temp_r"
  //    << temp_r << "\n";
  // }
  // End Debug Output

  // return{0.,0.,0.,0.};//debug

  return {2 * force_magnitude * displacement[0],
          2 * force_magnitude * displacement[1],
          2 * force_magnitude * displacement[2],
          0.0};  // 4th component is unused
}

InteractionForce* InteractionVelocity::NewCopy() const {
  return new InteractionVelocity();
}

}  // namespace bdm