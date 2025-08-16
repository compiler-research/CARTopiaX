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

#ifndef FORCES_TUMOR_CART_H_
#define FORCES_TUMOR_CART_H_
#include "core/interaction_force.h"
#include "core/operation/mechanical_forces_op.h"
#include "biodynamo.h"
#include "core/util/log.h"
#include "core/util/root.h"
#include "hyperparams.h"
#include "utils_aux.h"
#include "tumor_cell.h"

namespace bdm {

/**
 * @brief Custom interaction force implementation for velocity-based cell interactions
 * 
 * This class implements a specialized interaction force that takes into account
 * the velocity of cells when calculating forces between agents (tumor cells and CAR-T cells).
 * It extends the base InteractionForce class to provide custom force calculations
 * specific to the tumor-CAR-T cell interaction simulation.
 */
class InteractionVelocity : public InteractionForce {
 public:
  InteractionVelocity() = default;
  
  ~InteractionVelocity() override = default;

  /** @brief Calculate interaction force between two agents
   * 
   * Computes the force vector between two agents (cells) based on their
   * positions, properties, and velocities. This method is called by the
   * mechanical forces operation during each simulation step.
   * 
   * @param lhs Pointer to the first agent (left-hand side)
   * @param rhs Pointer to the second agent (right-hand side)
   * @return Real4 vector containing the force components (fx, fy, fz, magnitude)
   */
  Real4 Calculate(const Agent* lhs, const Agent* rhs) const override;

  InteractionForce* NewCopy() const override;
};

}  // namespace bdm

#endif  // FORCES_TUMOR_CART_H_