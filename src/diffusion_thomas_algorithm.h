// -----------------------------------------------------------------------------
// Copyright (C) 2025 Salvador de la Torre Gonzalez
// Co-author: Luciana Melina Luque
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------

#ifndef DIFFUSION_THOMAS_ALGORITHM_H_
#define DIFFUSION_THOMAS_ALGORITHM_H_

#include <utility>

#include "core/diffusion/diffusion_grid.h"

namespace bdm {

/** @brief Continuum model for the 3D heat equation with exponential decay
           \f$ \partial_t u = \nabla D \nabla u - \mu u \f$.
*/
class DiffusionThomasAlgorithm : public DiffusionGrid {
 public:
  DiffusionThomasAlgorithm() = default;
  DiffusionThomasAlgorithm(int substance_id, // Substance ID
                            std::string substance_name, // Substance name
                            real_t dc, // Diffusion coefficient
                            real_t mu,// Diffusion coefficient and decay constant
                            int resolution, // number of voxels in each direction
                            real_t dt,//time step
                            bool dirichlet_border);//if the border conditions are Dirichlet, this flag should not be necessary in a future version of BioDynaMo

  //Methonds to set new concentration at a given position
  void SetConcentration(real_t x, real_t y, real_t z, real_t amount){
    SetConcentration(GetBoxIndex(x, y, z), amount);
  };
  void SetConcentration(size_t idx, real_t amount);

  //These methods are overridden but empty because they are not used. This should be fixed in future versions of BioDynaMo
  void DiffuseWithClosedEdge(real_t dt) override{};
  void DiffuseWithOpenEdge(real_t dt) override{};
  void DiffuseWithNeumann(real_t dt) override{};
  void DiffuseWithPeriodic(real_t dt) override{};
  void DiffuseWithDirichlet(real_t dt) override{};


  //These methods are the important ones
  void DiffuseChemical(real_t dt);
  void Step(real_t dt) override;
  // ─────────────────────────────
  // Secretion or consumption of a substance following the differential equation
  // ∂ρ/∂t = ∇·(D ∇ρ) − λ · ρ + sum_{cells in voxel}((V_k / V_voxel) · [ S_k · ( ρ*_k − ρ ) − (S_k + U_k) · ρ ])
  // where:
  // ρ      = concentration of the substance in the microenvironment
  // S_k    = secretion rate of cell k
  // U_k    = uptake (consumption) rate of cell k
  // ρ*_k   = saturation (target) density for secretion
  // V_k    = volume of the cell k which is approximated to the default volume of a new tumor cell
  // V_voxel = volume of the voxel containing the cell
  // dt     = simulation time step
  //
  // In this class, we only model the secretion and consumption of the substance,
  // not its diffusion, which is:
  // (ρ − σ)/dt = sum_{cells in voxel}((V_k / V_voxel) · [ S_k · ( ρ*_k − ρ ) − (S_k + U_k) · ρ ])
  // where σ is the concentration of the substance in the voxel at the previous time step (it can include the diffusion term)
  // ρⁿ⁺¹ = (ρⁿ + dt · (V_k / V_voxel) · S_k · ρ*_k) 
  //        / [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  //
  // where:
  // ρⁿ     = current concentration
  // ρⁿ⁺¹   = updated concentration after secretion/uptake
  // This assumes secretion is toward a saturation level, and uptake is proportional to ρ
  // ─────────────────────────────

  void ComputeConsumptionsSecretions();
  
 private:
  size_t resolution_; // Number of voxels in each direction
  real_t d_space_; // Voxel side length in micrometers

  std::vector<real_t> thomas_denom_x_; // Denominators for x-direction
  std::vector<real_t> thomas_c_x_; // Coefficients for x-direction
  std::vector<real_t> thomas_denom_y_; // Denominators for y-direction
  std::vector<real_t> thomas_c_y_; // Coefficients for y-direction
  std::vector<real_t> thomas_denom_z_; // Denominators for z-direction
  std::vector<real_t> thomas_c_z_; // Coefficients for z-direction
  int jump_i_;
  int jump_j_;
  int jump_k_;

  real_t constant1_;
  real_t constant1a_;
  real_t constant2_;
  real_t constant3_;
  real_t constant3a_;

  bool dirichlet_border_; // Flag to indicate if the border conditions are Dirichlet

  // Function for initializing the coefficients of each direction
  void InitializeThomasAlgorithmVectors(std::vector<real_t>& thomas_denom,
                                        std::vector<real_t>& thomas_c);
  // Function to apply Dirichlet boundary conditions to the grid
  void ApplyDirichletBoundaryConditions();

  size_t GetBoxIndex(size_t x, size_t y, size_t z) const;
  BDM_CLASS_DEF_OVERRIDE(DiffusionThomasAlgorithm, 1);
};

}  // namespace bdm

#endif  // DIFFUSION_THOMAS_ALGORITHM_H_
