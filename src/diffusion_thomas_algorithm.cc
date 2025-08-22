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
#include "diffusion_thomas_algorithm.h"
#include "core/resource_manager.h"
#include "core/simulation.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "cart_cell.h"

namespace bdm {


DiffusionThomasAlgorithm::DiffusionThomasAlgorithm(int substance_id, std::string substance_name, real_t dc, real_t mu,int resolution, real_t dt, bool dirichlet_border)//time step
  : DiffusionGrid(substance_id, std::move(substance_name), dc, mu, resolution) {
  
  SetTimeStep(dt);
  //num of voxels in each direction
  resolution_ = GetResolution();
  // Voxel side length in micrometers
  d_space_ = kBoundedSpaceLength / resolution_; 

  dirichlet_border_ = dirichlet_border;

  jump_i_ = 1;
  jump_j_ = resolution_;
  jump_k_ = resolution_ * resolution_;

  //all diffusion coefficients are the same for all directions (isotropic)
  constant1_ = dc;
  constant1_ *=dt/(d_space_ * d_space_);
  constant1a_ = -constant1_;
  //decay constant
  constant2_ = mu;
  // Divide by 3 for the three directions
  constant2_ *= dt / 3.0; 

  constant3_ = 1.0 + 2 * constant1_ + constant2_;
  constant3a_ = 1.0 + constant1_ + constant2_;

  // Initialize the denominators and coefficients for the Thomas algorithm

  thomas_c_x_ = std::vector<real_t>(resolution_, constant1a_);
  thomas_denom_x_ = std::vector<real_t>(resolution_, constant3_);
  InitializeThomasAlgorithmVectors(thomas_denom_x_, thomas_c_x_);

  thomas_c_y_ = std::vector<real_t>(resolution_, constant1a_);
  thomas_denom_y_ = std::vector<real_t>(resolution_, constant3_);
  InitializeThomasAlgorithmVectors(thomas_denom_y_, thomas_c_y_);

  thomas_c_z_ = std::vector<real_t>(resolution_, constant1a_);
  thomas_denom_z_ = std::vector<real_t>(resolution_, constant3_);
  InitializeThomasAlgorithmVectors(thomas_denom_z_, thomas_c_z_);

}

void DiffusionThomasAlgorithm::InitializeThomasAlgorithmVectors(std::vector<real_t>& thomas_denom, std::vector<real_t>& thomas_c) {
  thomas_denom[0] = constant3a_;
  thomas_denom[resolution_ - 1] = constant3a_;
  if(resolution_ == 1) {
    thomas_denom[0] = 1.0 + constant2_;
  }
  thomas_c[0] /= thomas_denom[0];
  for (unsigned int i = 1; i < resolution_; ++i) {
    thomas_denom[i] += constant1_ * thomas_c[i - 1];
    thomas_c[i] /= thomas_denom[i];
  }
}

// Apply Dirichlet boundary conditions to the grid
void DiffusionThomasAlgorithm::ApplyDirichletBoundaryConditions() {
  real_t origin= GetDimensionsPtr()[0];
  real_t simulated_time = GetSimulatedTime();
  #pragma omp parallel
  {
    //We apply the Dirichlet boundary conditions to the first and last layers in each direction
    //For z=0 and z=resolution_-1
    #pragma omp for collapse(2)
    for (size_t y = 0; y < resolution_; y++) {
      for (size_t x = 0; x < resolution_; x++) {
        real_t real_x = origin + x * d_space_;
        real_t real_y = origin + y * d_space_;
        //For z=0
        size_t z=0;
        real_t real_z = origin + z * d_space_;
        SetConcentration(x, y, z, GetBoundaryCondition()->Evaluate(real_x, real_y, real_z, simulated_time));
        //For z=resolution_-1
        z = resolution_ - 1;
        real_z = origin + z * d_space_;
        SetConcentration(x, y, z, GetBoundaryCondition()->Evaluate(real_x, real_y, real_z, simulated_time));
      }
    }
    //For y=0 and y=resolution_-1
    #pragma omp for collapse(2)
    for (size_t z = 0; z < resolution_; z++) {
      for (size_t x = 0; x < resolution_; x++) {
        real_t real_x = origin + x * d_space_;
        real_t real_z = origin + z * d_space_;
        //For y=0
        size_t y=0;
        real_t real_y = origin + y * d_space_;
        SetConcentration(x, y, z, GetBoundaryCondition()->Evaluate(real_x, real_y, real_z, simulated_time));
        //For y=resolution_-1
        y = resolution_ - 1;
        real_y = origin + y * d_space_;
        SetConcentration(x, y, z, GetBoundaryCondition()->Evaluate(real_x, real_y, real_z, simulated_time));
      }
    }
    //For x=0 and x=resolution_-1
    #pragma omp for collapse(2)
    for (size_t z = 0; z < resolution_; z++) {
      for (size_t y = 0; y < resolution_; y++) {
        real_t real_y = origin + y * d_space_;
        real_t real_z = origin + z * d_space_;
        //For x=0
        size_t x=0;
        real_t real_x = origin + x * d_space_;
        SetConcentration(x, y, z, GetBoundaryCondition()->Evaluate(real_x, real_y, real_z, simulated_time));
        //For x=resolution_-1
        x = resolution_ - 1;
        real_x = origin + x * d_space_;
        SetConcentration(x, y, z, GetBoundaryCondition()->Evaluate(real_x, real_y, real_z, simulated_time));
      }
    }
  }

  
}

// Sets the concentration at a specific voxel
void DiffusionThomasAlgorithm::SetConcentration(size_t idx, real_t amount){
  ChangeConcentrationBy(idx, amount - GetAllConcentrations()[idx], InteractionMode::kAdditive, false);
};

// Flattens the 3D coordinates (x, y, z) into a 1D index
size_t DiffusionThomasAlgorithm::GetBoxIndex(size_t x, size_t y, size_t z) const {
  return z * resolution_ * resolution_ + y * resolution_ + x;
}

void DiffusionThomasAlgorithm::Step(real_t dt) {//instead of overwriting Step, in future versions of BioDynaMo, we should overwrite CheckParameters
  // check if diffusion coefficient and decay constant are 0
  // i.e. if we don't need to calculate diffusion update
  if (IsFixedSubstance()) {
    return;
  }
  DiffuseChemical(dt);

  //This should be done considering different border cases instead of using the dirichlet_border_ flag. However, there is a bug in BioDynaMo that makes bc_type be "Neumann" no matter what. In future versions of BioDynaMo this should be fixed

}

//This method solves the Diffusion Diferential equation using the Alternating Direction Implicit approach
void DiffusionThomasAlgorithm::DiffuseChemical(real_t dt) {


  // Change for the future: to add double buffer for paralelization

  if (dirichlet_border_) { ApplyDirichletBoundaryConditions();}

  //X-direction
  #pragma omp parallel for collapse(2)
  for( unsigned int k=0; k < resolution_; k++) {
    for( unsigned int j=0; j < resolution_; j++) {
      int ind = GetBoxIndex(0, j, k);

      SetConcentration(ind, GetAllConcentrations()[ind]/thomas_denom_x_[0]);
      // Forward elimination step for x direction
      for (unsigned int i = 1; i < resolution_ ; i++) {
        ind = GetBoxIndex(i, j, k);
        auto* all_concentrations = GetAllConcentrations();
        SetConcentration(ind, (all_concentrations[ind] + constant1_ * all_concentrations[ind-jump_i_]) / thomas_denom_x_[i]);
      }
      // Back substitution step for x direction
      for (int i = resolution_ - 2; i >= 0; i--) {
        ind = GetBoxIndex(i, j, k);
        auto* all_concentrations = GetAllConcentrations();
        SetConcentration(ind, all_concentrations[ind] - thomas_c_x_[i] * all_concentrations[ind + jump_i_]);
      }
    }
  }

  if (dirichlet_border_) { ApplyDirichletBoundaryConditions();}

  //Y-direction
  #pragma omp parallel for collapse(2)
  for( unsigned int k=0; k < resolution_; k++) {
    for( unsigned int i=0; i < resolution_; i++) {
      int ind = GetBoxIndex(i, 0, k);

      SetConcentration(ind, GetAllConcentrations()[ind]/thomas_denom_y_[0]);
      // Forward elimination step for y direction
      for (unsigned int j = 1; j < resolution_ ; j++) {
        ind = GetBoxIndex(i, j, k);
        auto* all_concentrations = GetAllConcentrations();
        SetConcentration(ind, (all_concentrations[ind] + constant1_ * all_concentrations[ind-jump_j_]) / thomas_denom_y_[j]);
      }
      // Back substitution step for y direction
      for (int j = resolution_ - 2; j >= 0; j--) {
        ind = GetBoxIndex(i, j, k);
        auto* all_concentrations = GetAllConcentrations();
        SetConcentration(ind, all_concentrations[ind] - thomas_c_y_[j] * all_concentrations[ind + jump_j_]);
      }
    }
  }

  if (dirichlet_border_) { ApplyDirichletBoundaryConditions();}

  //Z-direction
  #pragma omp parallel for collapse(2)
  for( unsigned int j=0; j < resolution_; j++) {
    for( unsigned int i=0; i < resolution_; i++) {
      int ind = GetBoxIndex(i, j, 0);
      SetConcentration(ind, GetAllConcentrations()[ind]/thomas_denom_z_[0]);
      // Forward elimination step for z direction
      for (unsigned int k = 1; k < resolution_ ; k++) {
        ind = GetBoxIndex(i, j, k);
        auto* all_concentrations = GetAllConcentrations();
        SetConcentration(ind, (all_concentrations[ind] + constant1_ * all_concentrations[ind-jump_k_]) / thomas_denom_z_[k]);
      }
      // Back substitution step for z direction
      for (int k = resolution_ - 2; k >= 0; k--) {
        ind = GetBoxIndex(i, j, k);
        auto* all_concentrations = GetAllConcentrations();
        SetConcentration(ind, all_concentrations[ind] - thomas_c_z_[k] * all_concentrations[ind + jump_k_]);
      }
    }
  }
  if (dirichlet_border_) { ApplyDirichletBoundaryConditions(); }

  //Change of concentration levels because of agents
  ComputeConsumptionsSecretions();

  return;


}


void DiffusionThomasAlgorithm::ComputeConsumptionsSecretions() {
  // This method is called to compute the consumptions and secretions of substances
  // by the tumor cells. It iterates over all agents and applies the consumption
  // and secretion behaviors defined in the TumorCell class.
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  real_t current_time = GetSimulatedTime();
  //in a future version of BioDynaMo this should be parallelized getting the agents inside each chemical voxel and trating each voxel independently.
  rm->ForEachAgent([this, current_time](bdm::Agent* agent) {
    if (auto* cell = dynamic_cast<TumorCell*>(agent)) {
      // Handle TumorCell agents
      const auto& pos = cell->GetPosition();
      real_t conc = this->GetValue(pos);
      real_t new_conc = cell->ConsumeSecreteSubstance(GetContinuumId(),conc);
      this->ChangeConcentrationBy(pos, new_conc - conc, InteractionMode::kAdditive, false);
    } else if (auto* cell = dynamic_cast<CartCell*>(agent)) {
      // Handle CartCell agents
      const auto& pos = cell->GetPosition();
      real_t conc = GetValue(pos);
      real_t new_conc = cell->ConsumeSecreteSubstance(GetContinuumId(),conc);
      ChangeConcentrationBy(pos, new_conc - conc, InteractionMode::kAdditive, false);
    }

  });

  return;
}

}  // namespace bdm
