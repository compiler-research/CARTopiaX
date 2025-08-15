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

#include "diffusion_thomas_algorithm.h"
#include "cart_cell.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "core/agent/agent.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace bdm {

DiffusionThomasAlgorithm::DiffusionThomasAlgorithm(int substance_id,
                                                   std::string substance_name,
                                                   real_t dc, real_t mu,
                                                   int resolution, real_t dt,
                                                   bool dirichlet_border)
    : DiffusionGrid(substance_id, std::move(substance_name), dc, mu,
                    resolution),
      resolution_(static_cast<int>(GetResolution())),
      d_space_(static_cast<real_t>(kBoundedSpaceLength) /
               static_cast<real_t>(resolution_)),
      dirichlet_border_(dirichlet_border),
      jump_i_(1),
      jump_j_(resolution_),
      jump_k_(resolution_ * resolution_),
      spatial_diffusion_coeff_(dc * dt / (d_space_ * d_space_)),
      neg_diffusion_factor_(-spatial_diffusion_coeff_),
      temporal_decay_coeff_(mu * dt / 3.0),
      central_coeff_(1.0 + 2 * spatial_diffusion_coeff_ +
                     temporal_decay_coeff_),
      edge_coeff_(1.0 + spatial_diffusion_coeff_ + temporal_decay_coeff_),
      thomas_c_x_(resolution_, neg_diffusion_factor_),
      thomas_denom_x_(resolution_, central_coeff_),
      thomas_c_y_(resolution_, neg_diffusion_factor_),
      thomas_denom_y_(resolution_, central_coeff_),
      thomas_c_z_(resolution_, neg_diffusion_factor_),
      thomas_denom_z_(resolution_, central_coeff_) {
  SetTimeStep(dt);

  // Initialize the denominators and coefficients for the Thomas algorithm
  InitializeThomasAlgorithmVectors(thomas_denom_x_, thomas_c_x_);
  InitializeThomasAlgorithmVectors(thomas_denom_y_, thomas_c_y_);
  InitializeThomasAlgorithmVectors(thomas_denom_z_, thomas_c_z_);
}

void DiffusionThomasAlgorithm::InitializeThomasAlgorithmVectors(
    std::vector<real_t>& thomas_denom, std::vector<real_t>& thomas_c) const {
  thomas_denom[0] = edge_coeff_;
  thomas_denom[resolution_ - 1] = edge_coeff_;
  if (resolution_ == 1) {
    thomas_denom[0] = 1.0 + temporal_decay_coeff_;
  }
  thomas_c[0] /= thomas_denom[0];
  for (int i = 1; i < resolution_; ++i) {
    thomas_denom[i] += spatial_diffusion_coeff_ * thomas_c[i - 1];
    thomas_c[i] /= thomas_denom[i];
  }
}

// Apply Dirichlet boundary conditions to the grid
void DiffusionThomasAlgorithm::ApplyDirichletBoundaryConditions() {
  // FIXME: Fix BioDynaMo by returning a view or c++20 std::span.
  const auto* dimensions_ptr = GetDimensionsPtr();
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const real_t origin = dimensions_ptr[0];
  const real_t simulated_time = GetSimulatedTime();
#pragma omp parallel
  {
// We apply the Dirichlet boundary conditions to the first and last layers in
// each direction For z=0 and z=resolution_-1
#pragma omp for collapse(2)
    for (int y = 0; y < resolution_; y++) {
      for (int x = 0; x < resolution_; x++) {
        const real_t real_x = origin + x * d_space_;
        const real_t real_y = origin + y * d_space_;
        // For z=0
        int z = 0;
        real_t real_z = origin + z * d_space_;
        SetConcentration(x, y, z,
                         GetBoundaryCondition()->Evaluate(
                             real_x, real_y, real_z, simulated_time));
        // For z=resolution_-1
        z = resolution_ - 1;
        real_z = origin + z * d_space_;
        SetConcentration(x, y, z,
                         GetBoundaryCondition()->Evaluate(
                             real_x, real_y, real_z, simulated_time));
      }
    }
// For y=0 and y=resolution_-1
#pragma omp for collapse(2)
    for (int z = 0; z < resolution_; z++) {
      for (int x = 0; x < resolution_; x++) {
        const real_t real_x = origin + x * d_space_;
        const real_t real_z = origin + z * d_space_;
        // For y=0
        int y = 0;
        real_t real_y = origin + y * d_space_;
        SetConcentration(x, y, z,
                         GetBoundaryCondition()->Evaluate(
                             real_x, real_y, real_z, simulated_time));
        // For y=resolution_-1
        y = resolution_ - 1;
        real_y = origin + y * d_space_;
        SetConcentration(x, y, z,
                         GetBoundaryCondition()->Evaluate(
                             real_x, real_y, real_z, simulated_time));
      }
    }
// For x=0 and x=resolution_-1
#pragma omp for collapse(2)
    for (int z = 0; z < resolution_; z++) {
      for (int y = 0; y < resolution_; y++) {
        const real_t real_y = origin + y * d_space_;
        const real_t real_z = origin + z * d_space_;
        // For x=0
        int x = 0;
        real_t real_x = origin + x * d_space_;
        SetConcentration(x, y, z,
                         GetBoundaryCondition()->Evaluate(
                             real_x, real_y, real_z, simulated_time));
        // For x=resolution_-1
        x = resolution_ - 1;
        real_x = origin + x * d_space_;
        SetConcentration(x, y, z,
                         GetBoundaryCondition()->Evaluate(
                             real_x, real_y, real_z, simulated_time));
      }
    }
  }
}

// Sets the concentration at a specific voxel
void DiffusionThomasAlgorithm::SetConcentration(size_t idx, real_t amount) {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const auto* all_concentrations = GetAllConcentrations();
  const real_t current_concentration = all_concentrations[idx];
  ChangeConcentrationBy(idx, amount - current_concentration,
                        InteractionMode::kAdditive, false);
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

// Flattens the 3D coordinates (x, y, z) into a 1D index
size_t DiffusionThomasAlgorithm::GetBoxIndex(size_t x, size_t y,
                                             size_t z) const {
  return z * resolution_ * resolution_ + y * resolution_ + x;
}

void DiffusionThomasAlgorithm::Step(real_t /*dt*/) {
  // check if diffusion coefficient and decay constant are 0
  // i.e. if we don't need to calculate diffusion update
  if (IsFixedSubstance()) {
    return;
  }
  DiffuseChemical();

  // This should be done considering different border cases instead of using the
  // dirichlet_border_ flag. However, there is a bug in BioDynaMo that makes
  // bc_type be "Neumann" no matter what. In future versions of BioDynaMo this
  // should be fixed
}

// This method solves the Diffusion Diferential equation using the Alternating
// Direction Implicit approach
void DiffusionThomasAlgorithm::DiffuseChemical() {
  ApplyBoundaryConditionsIfNeeded();

  // Solve for X-direction (direction = 0)
  SolveDirectionThomas(0);
  ApplyBoundaryConditionsIfNeeded();

  // Solve for Y-direction (direction = 1)
  SolveDirectionThomas(1);
  ApplyBoundaryConditionsIfNeeded();

  // Solve for Z-direction (direction = 2)
  SolveDirectionThomas(2);
  ApplyBoundaryConditionsIfNeeded();

  // Change of concentration levels because of agents
  ComputeConsumptionsSecretions();
}

void DiffusionThomasAlgorithm::ApplyBoundaryConditionsIfNeeded() {
  if (dirichlet_border_) {
    ApplyDirichletBoundaryConditions();
  }
}

void DiffusionThomasAlgorithm::SolveDirectionThomas(int direction) {
  const auto& thomas_denom = [this, direction]() -> const std::vector<real_t>& {
    if (direction == 0) {
      return thomas_denom_x_;
    }
    if (direction == 1) {
      return thomas_denom_y_;
    }
    // direction == 2
    return thomas_denom_z_;
  }();

  const auto& thomas_c = [this, direction]() -> const std::vector<real_t>& {
    if (direction == 0) {
      return thomas_c_x_;
    }
    if (direction == 1) {
      return thomas_c_y_;
    }
    // direction == 2
    return thomas_c_z_;
  }();

  const int jump = [this, direction]() -> int {
    if (direction == 0) {
      return jump_i_;
    }
    if (direction == 1) {
      return jump_j_;
    }
    // direction == 2
    return jump_k_;
  }();

#pragma omp parallel for collapse(2)
  for (int outer = 0; outer < resolution_; outer++) {
    for (int middle = 0; middle < resolution_; middle++) {
      // Forward elimination step
      ForwardElimination(direction, outer, middle, thomas_denom, jump);

      // Back substitution step
      BackSubstitution(direction, outer, middle, thomas_c, jump);
    }
  }
}

void DiffusionThomasAlgorithm::ForwardElimination(
    int direction, int outer, int middle,
    const std::vector<real_t>& thomas_denom, int jump) {
  // Get initial index based on direction
  size_t ind = GetLoopIndex(direction, outer, middle, 0);
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const auto* all_concentrations = GetAllConcentrations();
  const real_t initial_concentration = all_concentrations[ind];
  SetConcentration(ind, initial_concentration / thomas_denom[0]);

  // Forward elimination loop
  for (int inner = 1; inner < resolution_; inner++) {
    ind = GetLoopIndex(direction, outer, middle, inner);
    const real_t current_concentration = all_concentrations[ind];
    const real_t prev_concentration = all_concentrations[ind - jump];
    SetConcentration(ind, (current_concentration +
                           spatial_diffusion_coeff_ * prev_concentration) /
                              thomas_denom[inner]);
  }
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

void DiffusionThomasAlgorithm::BackSubstitution(
    int direction, int outer, int middle, const std::vector<real_t>& thomas_c,
    int jump) {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const auto* all_concentrations = GetAllConcentrations();

  // Back substitution loop
  for (int inner = resolution_ - 2; inner >= 0; inner--) {
    const size_t ind = GetLoopIndex(direction, outer, middle, inner);
    const real_t current_concentration = all_concentrations[ind];
    const real_t next_concentration = all_concentrations[ind + jump];
    SetConcentration(
        ind, current_concentration - thomas_c[inner] * next_concentration);
  }
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

size_t DiffusionThomasAlgorithm::GetLoopIndex(int direction, int outer,
                                              int middle, int inner) const {
  switch (direction) {
    case 0:  // X-direction: outer=k, middle=j, inner=i
      return GetBoxIndex(inner, middle, outer);
    case 1:  // Y-direction: outer=k, middle=i, inner=j
      return GetBoxIndex(middle, inner, outer);
    case 2:  // Z-direction: outer=j, middle=i, inner=k
      return GetBoxIndex(middle, outer, inner);
    default:
      return 0;
  }
}

void DiffusionThomasAlgorithm::ComputeConsumptionsSecretions() {
  // This method is called to compute the consumptions and secretions of
  // substances by the tumor cells. It iterates over all agents and applies the
  // consumption and secretion behaviors defined in the TumorCell class.
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // in a future version of BioDynaMo this should be parallelized getting the
  // agents inside each chemical voxel and treating each voxel independently.
  // Fixme:this dynamic casting could be optimized in a future version
  rm->ForEachAgent([this](bdm::Agent* agent) {
    if (auto* cell = dynamic_cast<TumorCell*>(agent)) {
      // Handle TumorCell agents
      const auto& pos = cell->GetPosition();
      const real_t conc = this->GetValue(pos);
      const real_t new_conc =
          cell->ConsumeSecreteSubstance(GetContinuumId(), conc);
      this->ChangeConcentrationBy(pos, new_conc - conc,
                                  InteractionMode::kAdditive, false);
    } else if (auto* cell = dynamic_cast<CarTCell*>(agent)) {
      // Handle CarTCell agents
      const auto& pos = cell->GetPosition();
      const real_t conc = GetValue(pos);
      const real_t new_conc =
          cell->ConsumeSecreteSubstance(GetContinuumId(), conc);
      ChangeConcentrationBy(pos, new_conc - conc, InteractionMode::kAdditive,
                            false);
    }
  });
}

}  // namespace bdm
