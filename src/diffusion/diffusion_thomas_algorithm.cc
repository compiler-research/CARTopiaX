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

#include "diffusion/diffusion_thomas_algorithm.h"
#include "interfaces/substance_interactor.h"
#include "params/hyperparams.h"
#include "core/agent/agent.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/param/param.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace bdm {

// Forward elimination for W contiguous lines at once (y/z sweeps).
// Line b's element i lives at c[base + i*jump + b]. Each line runs the
// same operations in the same order as the scalar version, so results
// are bit-identical; the W lines only overlap in execution.
template <int W>
inline void ForwardBlock(real_t* __restrict c, size_t base, int jump, int n,
                         const real_t* __restrict d, real_t k, real_t lo,
                         real_t hi) {
  real_t prev[W];
  for (int b = 0; b < W; ++b) {
    const real_t old = c[base + b];
    const real_t v = old / d[0];
    prev[b] = std::clamp(old + (v - old), lo, hi);
    c[base + b] = prev[b];
  }
  size_t ind = base;
  for (int i = 1; i < n; ++i) {
    ind += jump;
    for (int b = 0; b < W; ++b) {
      const real_t old = c[ind + b];
      const real_t v = (old + k * prev[b]) / d[i];
      prev[b] = std::clamp(old + (v - old), lo, hi);
      c[ind + b] = prev[b];
    }
  }
}
// Back substitution for W contiguous lines at once (y/z sweeps).
template <int W>
inline void BackBlock(real_t* __restrict c, size_t base, int jump, int n,
                      const real_t* __restrict tc, real_t lo, real_t hi) {
  real_t next[W];
  size_t ind = base + static_cast<size_t>(n - 1) * jump;
  for (int b = 0; b < W; ++b)
    next[b] = c[ind + b];  // last element
  for (int i = n - 2; i >= 0; --i) {
    ind -= jump;
    for (int b = 0; b < W; ++b) {
      const real_t old = c[ind + b];
      const real_t v = old - tc[i] * next[b];
      next[b] = std::clamp(old + (v - old), lo, hi);
      c[ind + b] = next[b];
    }
  }
}
// Forward elimination for W lines spaced 'ls' apart (x sweep).
// Neighbouring x lines are N voxels apart, so the values cannot be
// loaded contiguously; the win here is overlapping W independent
// dependency chains rather than SIMD.
template <int W>
inline void ForwardBlockStrided(real_t* __restrict c, size_t base, int jump,
                                int ls, int n, const real_t* __restrict d,
                                real_t k, real_t lo, real_t hi) {
  real_t prev[W];
  for (int b = 0; b < W; b++) {
    const size_t p = base + static_cast<size_t>(b) * ls;
    const real_t old = c[p];
    const real_t v = old / d[0];
    prev[b] = std::clamp(old + (v - old), lo, hi);
    c[p] = prev[b];
  }
  size_t ind = base;
  for (int i = 1; i < n; ++i) {
    ind += jump;
    for (int b = 0; b < W; b++) {
      const size_t p = ind + static_cast<size_t>(b) * ls;
      const real_t old = c[p];
      const real_t v = (old + k * prev[b]) / d[i];
      prev[b] = std::clamp(old + (v - old), lo, hi);
      c[p] = prev[b];
    }
  }
}
// Back substitution for W lines spaced 'ls' apart (x sweep).
template <int W>
inline void BackBlockStrided(real_t* __restrict c, size_t base, int jump,
                             int ls, int n, const real_t* __restrict tc,
                             real_t lo, real_t hi) {
  real_t next[W];
  size_t ind = base + static_cast<size_t>(n - 1) * jump;
  for (int b = 0; b < W; ++b)
    next[b] = c[ind + static_cast<size_t>(b) * ls];
  for (int i = n - 2; i >= 0; --i) {
    ind -= jump;
    for (int b = 0; b < W; ++b) {
      const size_t p = ind + static_cast<size_t>(b) * ls;
      const real_t old = c[p];
      const real_t v = old - tc[i] * next[b];
      next[b] = std::clamp(old + (v - old), lo, hi);
      c[p] = next[b];
    }
  }
}

DiffusionThomasAlgorithm::DiffusionThomasAlgorithm(int substance_id,
                                                   std::string substance_name,
                                                   real_t dc, real_t mu,
                                                   int resolution, real_t dt,
                                                   bool dirichlet_border)
    : DiffusionGrid(substance_id, std::move(substance_name), dc, mu,
                    resolution),
      resolution_(static_cast<int>(GetResolution())),
      d_space_(static_cast<real_t>(Simulation::GetActive()
                                       ->GetParam()
                                       ->Get<SimParam>()
                                       ->bounded_space_length) /
               static_cast<real_t>(resolution_)),
      dirichlet_border_(dirichlet_border),
      jump_i_(1),
      jump_j_(resolution_),
      jump_(resolution_ * resolution_),
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
// must be called from inside an omp parallel region; the loops
void DiffusionThomasAlgorithm::ApplyDirichletBoundaryConditions() {
  // FIXME: Fix BioDynaMo by returning a view or c++20 std::span.
  const int32_t* dimensions_ptr = GetDimensionsPtr();
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const real_t origin = dimensions_ptr[0];
  const real_t simulated_time = GetSimulatedTime();

// We apply the Dirichlet boundary conditions to the first and last layers in
// each direction For z=0 and z=resolution_-1
// THE USE OF nowait !!
// The three loops write different faces, so they need no barrier
// between them. The twelve cube edges are written by two loops each,
// with the same boundary value, so the overlap is benign.
#pragma omp for collapse(2)
  for (int y = 0; y < resolution_; y++) {
    for (int x = 0; x < resolution_; x++) {
      const real_t real_x = origin + x * d_space_;
      const real_t real_y = origin + y * d_space_;
      // For z=0
      int z = 0;
      real_t real_z = origin + z * d_space_;
      SetConcentration(x, y, z,
                       GetBoundaryCondition()->Evaluate(real_x, real_y, real_z,
                                                        simulated_time));
      // For z=resolution_-1
      z = resolution_ - 1;
      real_z = origin + z * d_space_;
      SetConcentration(x, y, z,
                       GetBoundaryCondition()->Evaluate(real_x, real_y, real_z,
                                                        simulated_time));
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
                       GetBoundaryCondition()->Evaluate(real_x, real_y, real_z,
                                                        simulated_time));
      // For y=resolution_-1
      y = resolution_ - 1;
      real_y = origin + y * d_space_;
      SetConcentration(x, y, z,
                       GetBoundaryCondition()->Evaluate(real_x, real_y, real_z,
                                                        simulated_time));
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

// Sets the concentration at a specific voxel
void DiffusionThomasAlgorithm::SetConcentration(size_t idx, real_t amount) {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const real_t* all_concentrations = GetAllConcentrations();
  const real_t current_concentration = all_concentrations[idx];
  ChangeConcentrationBy(idx, amount - current_concentration,
                        InteractionMode::kAdditive, false);
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

// Flattens the 3D coordinates (x, y, z) into a 1D index
size_t DiffusionThomasAlgorithm::GetBoxIndex(size_t x, size_t y,
                                             size_t z) const {
  assert(static_cast<int>(x) < resolution_ &&
         static_cast<int>(y) < resolution_ &&
         static_cast<int>(z) < resolution_ &&
         "GetBoxIndex: coordinate out of bounds");
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
// One parallel region spans all three sweeps and the boundary passes,
// instead of opening a region per sweep. The worksharing loops inside
// SolveDirectionThomas and ApplyDirichletBoundaryConditions bind to it.
// ComputeConsumptionsSecretions stays outside: it is serial and uses
// its own iteration over agents.
void DiffusionThomasAlgorithm::DiffuseChemical() {
#pragma omp parallel
  {
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
  }
  ComputeConsumptionsSecretions();
}

void DiffusionThomasAlgorithm::ApplyBoundaryConditionsIfNeeded() {
  if (dirichlet_border_) {
    ApplyDirichletBoundaryConditions();
  }
}

void DiffusionThomasAlgorithm::SolveDirectionThomas(int direction) {
  const std::array<const std::vector<real_t>*, 3> all_denoms = {
      &thomas_denom_x_, &thomas_denom_y_, &thomas_denom_z_};
  const std::array<const std::vector<real_t>*, 3> all_c = {
      &thomas_c_x_, &thomas_c_y_, &thomas_c_z_};
  const std::array<int, 3> all_jumps = {jump_i_, jump_j_, jump_};

  const std::vector<real_t>& thomas_denom = *all_denoms.at(direction);
  const std::vector<real_t>& thomas_c = *all_c.at(direction);
  const int jump = all_jumps.at(direction);
  const int n = resolution_;

  constexpr int kW = 8;
  const int tail_start = (n / kW) * kW;
  const int nblocks = tail_start / kW;
  real_t* c = const_cast<real_t*>(GetAllConcentrations());
  const real_t lo = GetLowerThreshold();
  const real_t hi = GetUpperThreshold();
  const real_t* d = thomas_denom.data();
  const real_t* tc = thomas_c.data();
  const real_t k = spatial_diffusion_coeff_;

  // must be called from inside an omp parallel region.
  if (direction == 0) {
    const int line_stride = jump_j_;  // = N

#pragma omp for collapse(2) nowait
    for (int outer = 0; outer < n; outer++) {
      for (int mblock = 0; mblock < nblocks; mblock++) {
        const int middle = mblock * kW;
        const size_t base = GetLoopIndex(direction, outer, middle, 0);
        ForwardBlockStrided<kW>(c, base, jump, line_stride, n, d, k, lo, hi);
        BackBlockStrided<kW>(c, base, jump, line_stride, n, tc, lo, hi);
      }
    }

#pragma omp for collapse(2)
    for (int outer = 0; outer < n; outer++) {
      for (int middle = tail_start; middle < n; middle++) {
        const size_t base = GetLoopIndex(direction, outer, middle, 0);
        ForwardBlockStrided<1>(c, base, jump, line_stride, n, d, k, lo, hi);
        BackBlockStrided<1>(c, base, jump, line_stride, n, tc, lo, hi);
      }
    }
  } else {
#pragma omp for collapse(2) nowait
    for (int outer = 0; outer < n; outer++) {
      for (int mblock = 0; mblock < nblocks; mblock++) {
        const int middle = mblock * kW;
        const size_t base = GetLoopIndex(direction, outer, middle, 0);
        ForwardBlock<kW>(c, base, jump, n, d, k, lo, hi);
        BackBlock<kW>(c, base, jump, n, tc, lo, hi);
      }
    }

#pragma omp for collapse(2)
    for (int outer = 0; outer < n; outer++) {
      for (int middle = tail_start; middle < n; middle++) {
        const size_t base = GetLoopIndex(direction, outer, middle, 0);
        ForwardBlock<1>(c, base, jump, n, d, k, lo, hi);
        BackBlock<1>(c, base, jump, n, tc, lo, hi);
      }
    }
  }
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
  ResourceManager* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // in a future version of BioDynaMo this should be parallelized getting the
  // agents inside each chemical voxel and treating each voxel independently.

  rm->ForEachAgent([this](bdm::Agent* agent) {
    if (auto* interactor = dynamic_cast<ISubstanceInteractor*>(agent)) {
      const Real3& pos = agent->GetPosition();
      const real_t conc = GetValue(pos);
      const real_t new_conc =
          interactor->ConsumeSecreteSubstance(GetContinuumId(), conc);
      ChangeConcentrationBy(pos, new_conc - conc, InteractionMode::kAdditive,
                            false);
    }
  });
}

}  // namespace bdm
