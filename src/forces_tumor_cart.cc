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

#include "forces_tumor_cart.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "utils_aux.h"
#include "core/agent/agent.h"
#include "core/agent/cell.h"
#include "core/container/math_array.h"
#include "core/interaction_force.h"
#include "core/real_t.h"
#include <algorithm>
#include <cmath>
#include <memory>

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

  const real_t radius_a = a->GetDiameter() / kHalf;
  const real_t radius_b = b->GetDiameter() / kHalf;
  const real_t combined_radius = radius_a + radius_b;
  real_t temp_r = 0.0;

  const auto* a_tumor = dynamic_cast<const TumorCell*>(a);
  const auto* b_tumor = dynamic_cast<const TumorCell*>(b);

  if (distance < combined_radius) {
    // 1 - d/combined_radius
    temp_r = 1.0 - distance / combined_radius;
    // (1 - d/combined_radius)^2
    temp_r *= temp_r;

    real_t repulsion = NAN;

    if ((a_tumor != nullptr) && (b_tumor != nullptr)) {
      // two tumor cells
      // std::sqrt(kRepulsionTumorTumor * kRepulsionTumorTumor);
      repulsion = kRepulsionTumorTumor;
    } else if ((a_tumor == nullptr) && (b_tumor == nullptr)) {
      // two CAR-T cells
      // std::sqrt(kRepulsionCartCart*kRepulsionCartCart);
      repulsion = kRepulsionCartCart;
    } else {
      // one tumor cell and one CAR-T
      repulsion = std::sqrt(kRepulsionCartTumor * kRepulsionTumorCart);
    }

    temp_r *= repulsion;
  }

  // Adhesion
  const real_t max_interaction_distance =
      kMaxRelativeAdhesionDistance * combined_radius;

  if (distance < max_interaction_distance) {
    // 1 - d/S
    real_t temp_a = 1.0 - distance / max_interaction_distance;
    // (1-d/S)^2
    temp_a *= temp_a;
    // Initialize to NAN
    real_t adhesion = NAN;
    if ((a_tumor != nullptr) && (b_tumor != nullptr)) {
      // two tumor cells
      adhesion = kAdhesionTumorTumor;
    } else if ((a_tumor == nullptr) && (b_tumor == nullptr)) {
      // two CAR-T cells
      adhesion = kAdhesionCartCart;
    } else {
      // one tumor cell and one CAR-T
      adhesion = std::sqrt(kAdhesionCartTumor * kAdhesionTumorCart);
    }

    temp_a *= adhesion;
    temp_r -= temp_a;
  }

  if (std::abs(temp_r) < kEpsilon) {
    return {0.0, 0.0, 0.0, 0.0};
  }
  const real_t force_magnitude = temp_r / distance;

  return {force_magnitude * displacement[0], force_magnitude * displacement[1],
          force_magnitude * displacement[2],
          // 4th component is unused
          0.0};
}

InteractionForce* InteractionVelocity::NewCopy() const {
  return std::make_unique<InteractionVelocity>().release();
}

}  // namespace bdm