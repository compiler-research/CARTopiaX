/*
 * Copyright 2026 compiler-research.org, Salvador de la Torre Gonzalez
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

#include "diffusion/cylinder_wall_boundary_condition.h"

namespace bdm {

CylinderWallBoundaryCondition::CylinderWallBoundaryCondition(real_t value,
                                                              real_t min_z,
                                                              real_t max_z)
    : value_(value), min_z_(min_z), max_z_(max_z) {}

real_t CylinderWallBoundaryCondition::Evaluate(real_t /*x*/, real_t /*y*/,
                                               real_t z,
                                               real_t /*time*/) const {
  if (z < min_z_ || z > max_z_) {
    return 0.0;
  }
  return value_;
}

}  // namespace bdm