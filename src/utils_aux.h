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


#ifndef CORE_UTIL_UTILS_AUX_H_
#define CORE_UTIL_UTILS_AUX_H_

#include "biodynamo.h"
#include "core/util/random.h"
#include "hyperparams.h"
#include "tumor_cell.h"

namespace bdm {
    class TumorCell;  // Forward declaration

// Samples a Gaussian value with given mean and standard deviation but all negative values are mapped to zero
real_t SamplePositiveGaussian(float mean, float sigma);

// Samples a random antigen pattern from the predefined antigen patterns for CAR-T cells
// inline std::map<std::string, bool> SampleAntigenPattern(const std::vector<std::pair<std::map<std::string, bool>, float>>& possibilitiesRecognizedAntigensCart) {
//     float accum = 0.0f;

//     float rnumber = Simulation::GetActive()->GetRandom()->Uniform(0.0f, 1.0f);

//     for (const auto& [dictionary, probability] : possibilitiesRecognizedAntigensCart) {
//         accum += probability;
//         if (rnumber <= accum) {
//             return dictionary;
//         }
//     }

//     // Fallback in case no pattern is selected (should not happen)
//     return possibilitiesRecognizedAntigensCart.back().first;
// }

// ─────────────────────────────
// Behavior: Secretion or consumption of a substance following the differential equation
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

// class ConsumptionSecretion : public Behavior {
//   BDM_BEHAVIOR_HEADER(ConsumptionSecretion, Behavior, 1);

//  public:
//   ConsumptionSecretion() = default;
//   explicit ConsumptionSecretion(const std::string& substance, real_t quantity_consumption, real_t quantity_secretion, real_t substance_saturation_density);
//   explicit ConsumptionSecretion(DiffusionGrid* dgrid, real_t quantity_consumption, real_t quantity_secretion,real_t substance_saturation_density);

//   virtual ~ConsumptionSecretion() = default;

//   void Initialize(const NewAgentEvent& event) override;

//   void Run(Agent* agent) override;

//   real_t GetConsumption() const { return quantity_consumption_; }
//   real_t GetSecretion() const { return quantity_secretion_; }

//   void SetQuantities(real_t quantity_consumption, real_t quantity_secretion);// Set quantities for consumption and secretion by giving them already scaled

//   void ComputeConstants(real_t total_volume);

//   private:
//     DiffusionGrid* dgrid_ = nullptr;
//     real_t constant1_;//  = dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation  =  dt · (V_k / V_voxel) · S_k · ρ*_k) 
//     real_t constant2_;// = 1 + dt*(cell_volume/voxel_volume)*(quantity_secretion + quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
//     real_t quantity_consumption_;
//     real_t quantity_secretion_;
//     real_t substance_saturation_density_; // Saturation density for the substance, the agent tries to secrete till this value
// };

std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius);

//Function to compute the number of tumor cells of each type and the radius of the tumor
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, real_t> ComputeNumberTumorCellsAndRadius();

// Function to output summary CSV
struct OutputSummary : public StandaloneOperationImpl {
  BDM_OP_HEADER(OutputSummary);

  uint64_t frequency_ = 1;

  void operator()() override;
};

// Register with CPU as compute target
inline BDM_REGISTER_OP(OutputSummary, "OutputSummary", kCpu);

}  // namespace bdm

#endif  // CORE_UTIL_UTILS_AUX_H_