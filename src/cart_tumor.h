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

#ifndef CART_TUMOR_H_
#define CART_TUMOR_H_

#include "biodynamo.h"
#include "tumor_cell.h"
#include "cart_cell.h"
#include "diffusion_thomas_algorithm.h"
#include "forces_tumor_cart.h" 
#include "core/environment/uniform_grid_environment.h"
#include "core/operation/mechanical_forces_op.h"
namespace bdm {



// List the extracellular substances
enum Substances { kImmunostimulatoryFactor, kOxygen };

inline int Simulate(int argc, const char** argv) {



  // Set simulation bounds
  auto set_param = [](Param* param) {
    param->random_seed = kSeed; // Set a fixed random seed for reproducibility
    param->bound_space = Param::BoundSpaceMode::kTorus;// Periodic boundary
    param->min_bound = -kBoundedSpaceLength / 2;
    param->max_bound = kBoundedSpaceLength/2;  // Cube of 1000x1000x1000 centered at origin
    param->simulation_time_step = kDt;
  };
  


  Simulation simulation(argc, argv, set_param);
  auto* ctxt = simulation.GetExecutionContext();

  //Change Forces
  auto* scheduler = simulation.GetScheduler();

  auto* op = scheduler->GetOps("mechanical forces")[0];
  op->GetImplementation<MechanicalForcesOp>()->SetInteractionForce(new InteractionVelocity());

  auto* env = dynamic_cast<UniformGridEnvironment*>(Simulation::GetActive()->GetEnvironment());
  env->SetBoxLength(kLengthBoxMechanics); // Fix the box length for the uniform grid environment

  // ───────────────────────────────────────
  // Define Substances
  // ───────────────────────────────────────
  auto* rm = Simulation::GetActive()->GetResourceManager();

  // Oxygen
  // substance_id, name, diffusion_coefficient, decay_constant, resolution, time_step
  auto* oxygen_grid = new DiffusionThomasAlgorithm(
      kOxygen, "oxygen",
      kDiffusionCoefficientOxygen,// 100000 micrometers^2/minute
      kDecayConstantOxygen, // 0.1 minutes^-1
      kResolutionGridSubstances,
      kDtSubstances,
      true); // true indicates Dirichlet border conditions
  rm->AddContinuum(oxygen_grid);

  // Immunostimulatory Factor
  // substance_id, name, diffusion_coefficient, decay_constant, resolution
  auto* immunostimulatory_factor_grid = new DiffusionThomasAlgorithm(
      kImmunostimulatoryFactor, "immunostimulatory_factor",
      kDiffusionCoefficientImmunostimulatoryFactor, // 1000 micrometers^2/minute
      kDecayConstantImmunostimulatoryFactor, // 0.016 minutes^-1
      kResolutionGridSubstances,
      kDtSubstances,
    false); // false indicates Neumann border conditions
  rm->AddContinuum(immunostimulatory_factor_grid);


  
  // Boundary Conditions Dirichlet: simulating absorption or total loss at the boundaries of the space.
  //Oxygen comming from the borders (capillary vessels)
  ModelInitializer::AddBoundaryConditions(
    kOxygen, BoundaryConditionType::kDirichlet,
    std::make_unique<ConstantBoundaryCondition>(kOxygenReferenceLevel));// kOxygenReferenceLevel mmHg is the physiological level of oxygen in tissues, o2 saturation is 100% at this level

  //This is useless now but should be added this way in a future version of BioDynaMo
  ModelInitializer::AddBoundaryConditions(
      kImmunostimulatoryFactor, BoundaryConditionType::kNeumann, nullptr);

  //Initialize oxygen voxels
  ModelInitializer::InitializeSubstance(kOxygen, [](real_t x, real_t y, real_t z) {
    return kInitialOxygenLevel; // Set all voxels to kInitialOxygenLevel mmHg
  });

  // ───────────────────────────────────────
  // One spherical tumor of radius kInitialRadiusTumor in the center of the simulation space
  // ───────────────────────────────────────
  std::vector<Real3> positions=CreateSphereOfTumorCells(kInitialRadiusTumor);//positions of the cells
  // positions={{0.1,0.1,0.1},{19.9,19.9,19.9}};//Debug {0.1,0.1,0.1},{19.9,19.9,19.9}
  for (const auto& pos : positions) {
    TumorCell* tumor_cell = new TumorCell(pos);
    tumor_cell->AddBehavior(new StateControlGrowProliferate());
    ctxt->AddAgent(tumor_cell);
  }


  // //debug
  // CartCell* cart_cell = new CartCell({0.,0.,0.});
  // cart_cell->AddBehavior(new StateControlCart());
  // ctxt->AddAgent(cart_cell);




  //OutputSummary operation
  auto* summary_op = new bdm::Operation("OutputSummary");
  summary_op->frequency_ = kOutputCsvInterval; // Set the interval for outputting CSV files
  summary_op->AddOperationImpl(bdm::kCpu, new bdm::OutputSummary());
  scheduler->ScheduleOp(summary_op);



  // ───────────────────────────────────────
  // Run simulation
  // ───────────────────────────────────────
  scheduler->Simulate(1+kTotalMinutesToSimulate/kDt);//simulate kTotalMinutesToSimulate minutes including the last minute 
  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // CART_TUMOR_H_
