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

#ifndef CYLINDER_WALL_BOUNDARY_CONDITION_H_
#define CYLINDER_WALL_BOUNDARY_CONDITION_H_

#include "core/diffusion/diffusion_grid.h"
#include "core/real_t.h"
#include "core/util/root.h"

namespace bdm {

/// @brief Dirichlet boundary condition that secretes a substance wherever
/// z falls within [min_z_, max_z_], regardless of which face of the domain
/// is being evaluated (x, y, or z faces). Outside that range, the boundary
/// condition evaluates to 0.
///
/// In the typical use case — modeling secretion through the side walls of
/// a cylindrical capillary region — min_z_ and max_z_ are set strictly
/// inside the domain, so the top and bottom faces (z = domain min / domain
/// max) naturally fall outside [min_z_, max_z_] and evaluate to 0.
///
/// However, this is a consequence of the chosen z-range, not a hardcoded
/// exclusion: if min_z_ is set exactly to the domain's minimum z (the
/// floor), the floor will also secrete value_ at that layer. Likewise, if
/// max_z_ is set exactly to the domain's maximum z (the roof), the roof
/// will also secrete value_ at that layer. Callers who need to guarantee
/// that floor/roof never secrete, regardless of the configured range,
/// should choose min_z_/max_z_ strictly inside the domain bounds.
class CylinderWallBoundaryCondition : public BoundaryCondition {
  using BoundaryCondition::BoundaryCondition;

 public:
  /// @param value Constant value secreted when min_z <= z <= max_z
  /// @param min_z Lower bound (inclusive) of the z-range where value_ is
  /// secreted. If set to the domain's minimum z, the floor also secretes.
  /// @param max_z Upper bound (inclusive) of the z-range where value_ is
  /// secreted. If set to the domain's maximum z, the roof also secretes.
  CylinderWallBoundaryCondition(real_t value, real_t min_z, real_t max_z);

  /// @brief see BoundaryCondition::Evaluate()
  real_t Evaluate(real_t x, real_t y, real_t z, real_t time) const final;

 private:
  /// Constant value of the boundary condition within [min_z_, max_z_]
  real_t value_ = 0.0;
  /// Lower z-bound (inclusive) where the boundary condition is active
  real_t min_z_ = 0.0;
  /// Upper z-bound (inclusive) where the boundary condition is active
  real_t max_z_ = 0.0;

  BDM_CLASS_DEF_OVERRIDE(CylinderWallBoundaryCondition, 1);
};

}  // namespace bdm

#endif  // CORE_DIFFUSION_CYLINDER_WALL_BOUNDARY_CONDITION_H_
