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

#ifndef DIFFUSION_THOMAS_ALGORITHM_H_
#define DIFFUSION_THOMAS_ALGORITHM_H_

#include "core/diffusion/diffusion_grid.h"
#include "core/real_t.h"
#include "core/util/root.h"
#include <cstddef>
#include <string>
#include <vector>

namespace bdm {

/// Continuum model for the 3D heat equation with exponential decay
///
/// Implements the diffusion equation, solved implicitly: ∂t u = ∇D∇u - μu
/// Uses the Thomas algorithm for solving tridiagonal systems efficiently.
class DiffusionThomasAlgorithm : public DiffusionGrid {
 public:
  DiffusionThomasAlgorithm()
      : resolution_(0),
        d_space_(0.0),
        dirichlet_border_(false),
        jump_i_(0),
        jump_j_(0),
        jump_k_(0),
        spatial_diffusion_coeff_(0.0),
        neg_diffusion_factor_(0.0),
        temporal_decay_coeff_(0.0),
        central_coeff_(0.0),
        edge_coeff_(0.0) {}

  DiffusionThomasAlgorithm(int substance_id, std::string substance_name,
                           real_t dc, real_t mu, int resolution, real_t dt,
                           bool dirichlet_border);

  /// Concentration setters
  void SetConcentration(int x, int y, int z, real_t amount) {
    SetConcentration(GetBoxIndex(x, y, z), amount);
  }

  void SetConcentration(size_t idx, real_t amount);

  /// These methods are overridden but empty because they are not used.
  /// This should be fixed in future versions of BioDynaMo.
  void DiffuseWithClosedEdge(real_t /*dt*/) override {}
  void DiffuseWithOpenEdge(real_t /*dt*/) override {}
  void DiffuseWithNeumann(real_t /*dt*/) override {}
  void DiffuseWithPeriodic(real_t /*dt*/) override {}
  void DiffuseWithDirichlet(real_t /*dt*/) override {}

  /// Perform chemical diffusion using Thomas algorithm
  ///
  /// Computes the diffusion of the substance using the Thomas algorithm
  /// for solving tridiagonal systems efficiently.
  void DiffuseChemical();

  /// Execute one simulation step
  ///
  /// Main stepping function that performs one time step of the simulation,
  /// including diffusion and cellular consumption/secretion.
  ///
  /// @param dt Time step for the simulation
  void Step(real_t dt) override;

  /// Compute cellular consumption and secretion effects
  ///
  /// Handles secretion or consumption of substances following the differential
  /// equation:
  ///
  /// ∂ρ/∂t = ∇·(D ∇ρ) − λ · ρ + sum_{cells in voxel}((V_k / V_voxel) · [ S_k ·
  /// ( ρ*_k − ρ ) − (S_k + U_k) · ρ ])
  ///
  /// Where:
  /// - ρ      = concentration of the substance in the microenvironment
  /// - S_k    = secretion rate of cell k
  /// - U_k    = uptake (consumption) rate of cell k
  /// - ρ*_k   = saturation (target) density for secretion
  /// - V_k    = volume of cell k (approximated to default volume of new tumor
  /// cell)
  /// - V_voxel = volume of the voxel containing the cell
  /// - dt     = simulation time step
  ///
  /// In this class, we only model the secretion and consumption of the
  /// substance, not its diffusion, which follows: (ρ − σ)/dt = sum_{cells in
  /// voxel}((V_k / V_voxel) · [ S_k · ( ρ*_k − ρ ) − (S_k + U_k) · ρ ])
  ///
  /// Where σ is the concentration at the previous time step (may include
  /// diffusion term). The solution is: ρⁿ⁺¹ = (ρⁿ + dt · (V_k / V_voxel) · S_k
  /// · ρ*_k) / [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  ///
  /// Where:
  /// - ρⁿ     = current concentration
  /// - ρⁿ⁺¹   = updated concentration after secretion/uptake
  ///
  /// This assumes secretion is toward a saturation level, and uptake is
  /// proportional to ρ.
  ///
  /// In a future version, consider using a Behavior associated to each agent
  /// but controlling the time in which it is applied so that it is executed
  /// always after the diffusion module
  void ComputeConsumptionsSecretions();

 private:
  /// Number of voxels in each spatial direction
  int resolution_;

  /// Voxel side length in micrometers
  real_t d_space_;

  /// Flag indicating Dirichlet boundary conditions
  bool dirichlet_border_;

  /// Index jump for i-direction (x-axis)
  int jump_i_;

  /// Index jump for j-direction (y-axis)
  int jump_j_;

  /// Index jump for k-direction (z-axis)
  int jump_k_;

  /// First diffusion constant
  real_t spatial_diffusion_coeff_;

  /// Alternative first diffusion constant
  real_t neg_diffusion_factor_;

  /// Second diffusion constant
  real_t temporal_decay_coeff_;

  /// Third diffusion constant
  real_t central_coeff_;

  /// Alternative third diffusion constant
  real_t edge_coeff_;

  /// Coefficients for x-direction Thomas algorithm
  std::vector<real_t> thomas_c_x_;

  /// Denominators for x-direction Thomas algorithm
  std::vector<real_t> thomas_denom_x_;

  /// Coefficients for y-direction Thomas algorithm
  std::vector<real_t> thomas_c_y_;

  /// Denominators for y-direction Thomas algorithm
  std::vector<real_t> thomas_denom_y_;

  /// Coefficients for z-direction Thomas algorithm
  std::vector<real_t> thomas_c_z_;

  /// Denominators for z-direction Thomas algorithm
  std::vector<real_t> thomas_denom_z_;

  /// Initialize Thomas algorithm coefficient vectors
  ///
  /// Sets up the precomputed coefficients for efficient Thomas algorithm
  /// execution in the specified direction.
  ///
  /// @param thomas_denom Reference to denominator vector to initialize
  /// @param thomas_c Reference to coefficient vector to initialize
  void InitializeThomasAlgorithmVectors(std::vector<real_t>& thomas_denom,
                                        std::vector<real_t>& thomas_c) const;

  /// Apply Dirichlet boundary conditions to the diffusion grid
  ///
  /// Sets the boundary values according to Dirichlet boundary conditions,
  /// maintaining constant values at the grid boundaries.
  void ApplyDirichletBoundaryConditions();

  /// Convert 3D coordinates to linear index
  ///
  /// @param x X-coordinate in voxel space
  /// @param y Y-coordinate in voxel space
  /// @param z Z-coordinate in voxel space
  /// @return Linear index in the flattened 3D array
  size_t GetBoxIndex(size_t x, size_t y, size_t z) const;

  /// Apply boundary conditions if Dirichlet boundaries are enabled
  void ApplyBoundaryConditionsIfNeeded();

  /// Solve Thomas algorithm for a specific direction
  ///
  /// @param direction Direction to solve (0=X, 1=Y, 2=Z)
  void SolveDirectionThomas(int direction);

  /// Perform forward elimination step of Thomas algorithm
  ///
  /// @param direction Direction being solved (0=X, 1=Y, 2=Z)
  /// @param outer Outer loop index
  /// @param middle Middle loop index
  /// @param thomas_denom Precomputed denominators for this direction
  /// @param jump Index jump value for this direction
  void ForwardElimination(int direction, int outer, int middle,
                          const std::vector<real_t>& thomas_denom, int jump);

  /// Perform back substitution step of Thomas algorithm
  ///
  /// @param direction Direction being solved (0=X, 1=Y, 2=Z)
  /// @param outer Outer loop index
  /// @param middle Middle loop index
  /// @param thomas_c Precomputed coefficients for this direction
  /// @param jump Index jump value for this direction
  void BackSubstitution(int direction, int outer, int middle,
                        const std::vector<real_t>& thomas_c, int jump);

  /// Get the linear index for given direction and loop indices
  ///
  /// @param direction Direction (0=X, 1=Y, 2=Z)
  /// @param outer Outer loop index
  /// @param middle Middle loop index
  /// @param inner Inner loop index
  /// @return Linear index in the flattened array
  size_t GetLoopIndex(int direction, int outer, int middle, int inner) const;

  BDM_CLASS_DEF_OVERRIDE(DiffusionThomasAlgorithm, 1);
};

}  // namespace bdm

#endif  // DIFFUSION_THOMAS_ALGORITHM_H_
