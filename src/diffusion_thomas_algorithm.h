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

#ifndef DIFFUSION_THOMAS_ALGORITHM_H_
#define DIFFUSION_THOMAS_ALGORITHM_H_

#include <utility>

#include "core/diffusion/diffusion_grid.h"

namespace bdm {

/// Continuum model for the 3D heat equation with exponential decay
/// 
/// Implements the diffusion equation, solved implicitly: ∂t u = ∇D∇u - μu
/// Uses the Thomas algorithm for solving tridiagonal systems efficiently.
class DiffusionThomasAlgorithm : public DiffusionGrid {
 public:
  DiffusionThomasAlgorithm() = default;
  
  DiffusionThomasAlgorithm(int substance_id,
                            std::string substance_name,
                            real_t dc,
                            real_t mu,
                            int resolution,
                            real_t dt,
                            bool dirichlet_border);

  /// Concentration setters
  void SetConcentration(real_t x, real_t y, real_t z, real_t amount){
    SetConcentration(GetBoxIndex(x, y, z), amount);
  };
  
  void SetConcentration(size_t idx, real_t amount);


  /// These methods are overridden but empty because they are not used.
  /// This should be fixed in future versions of BioDynaMo.
  void DiffuseWithClosedEdge(real_t dt) override{};
  void DiffuseWithOpenEdge(real_t dt) override{};
  void DiffuseWithNeumann(real_t dt) override{};
  void DiffuseWithPeriodic(real_t dt) override{};
  void DiffuseWithDirichlet(real_t dt) override{};


  /// Perform chemical diffusion using Thomas algorithm
  /// 
  /// Computes the diffusion of the substance using the Thomas algorithm
  /// for solving tridiagonal systems efficiently.
  /// 
  /// @param dt Time step for the diffusion computation
  void DiffuseChemical(real_t dt);
  
  /// Execute one simulation step
  /// 
  /// Main stepping function that performs one time step of the simulation,
  /// including diffusion and cellular consumption/secretion.
  /// 
  /// @param dt Time step for the simulation
  void Step(real_t dt) override;

  /// Compute cellular consumption and secretion effects
  /// 
  /// Handles secretion or consumption of substances following the differential equation:
  /// 
  /// ∂ρ/∂t = ∇·(D ∇ρ) − λ · ρ + sum_{cells in voxel}((V_k / V_voxel) · [ S_k · ( ρ*_k − ρ ) − (S_k + U_k) · ρ ])
  /// 
  /// Where:
  /// - ρ      = concentration of the substance in the microenvironment
  /// - S_k    = secretion rate of cell k
  /// - U_k    = uptake (consumption) rate of cell k
  /// - ρ*_k   = saturation (target) density for secretion
  /// - V_k    = volume of cell k (approximated to default volume of new tumor cell)
  /// - V_voxel = volume of the voxel containing the cell
  /// - dt     = simulation time step
  /// 
  /// In this class, we only model the secretion and consumption of the substance,
  /// not its diffusion, which follows:
  /// (ρ − σ)/dt = sum_{cells in voxel}((V_k / V_voxel) · [ S_k · ( ρ*_k − ρ ) − (S_k + U_k) · ρ ])
  /// 
  /// Where σ is the concentration at the previous time step (may include diffusion term).
  /// The solution is:
  /// ρⁿ⁺¹ = (ρⁿ + dt · (V_k / V_voxel) · S_k · ρ*_k) / [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  /// 
  /// Where:
  /// - ρⁿ     = current concentration
  /// - ρⁿ⁺¹   = updated concentration after secretion/uptake
  /// 
  /// This assumes secretion is toward a saturation level, and uptake is proportional to ρ.
  /// 
  /// In a future version, consider using a Behavior associated to each agent but controlling the time in which it is applied so that it is executed always after the diffusion module
  void ComputeConsumptionsSecretions();

 private:
  /// Number of voxels in each spatial direction
  size_t resolution_;
  
  /// Voxel side length in micrometers
  real_t d_space_;

  /// Denominators for x-direction Thomas algorithm
  std::vector<real_t> thomas_denom_x_;
  
  /// Coefficients for x-direction Thomas algorithm
  std::vector<real_t> thomas_c_x_;
  
  /// Denominators for y-direction Thomas algorithm
  std::vector<real_t> thomas_denom_y_;
  
  /// Coefficients for y-direction Thomas algorithm
  std::vector<real_t> thomas_c_y_;
  
  /// Denominators for z-direction Thomas algorithm
  std::vector<real_t> thomas_denom_z_;
  
  /// Coefficients for z-direction Thomas algorithm
  std::vector<real_t> thomas_c_z_;
  
  /// Index jump for i-direction (x-axis)
  int jump_i_;
  
  /// Index jump for j-direction (y-axis)
  int jump_j_;
  
  /// Index jump for k-direction (z-axis)
  int jump_k_;
  
  /// First diffusion constant
  real_t constant1_;
  
  /// Alternative first diffusion constant
  real_t constant1a_;
  
  /// Second diffusion constant
  real_t constant2_;
  
  /// Third diffusion constant
  real_t constant3_;
  
  /// Alternative third diffusion constant
  real_t constant3a_;
  
  /// Flag indicating Dirichlet boundary conditions
  bool dirichlet_border_;

  /// Initialize Thomas algorithm coefficient vectors
  /// 
  /// Sets up the precomputed coefficients for efficient Thomas algorithm
  /// execution in the specified direction.
  /// 
  /// @param thomas_denom Reference to denominator vector to initialize
  /// @param thomas_c Reference to coefficient vector to initialize
  void InitializeThomasAlgorithmVectors(std::vector<real_t>& thomas_denom,
                                        std::vector<real_t>& thomas_c);
  
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


  BDM_CLASS_DEF_OVERRIDE(DiffusionThomasAlgorithm, 1);
};

}  // namespace bdm

#endif  // DIFFUSION_THOMAS_ALGORITHM_H_
