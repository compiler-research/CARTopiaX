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

#include <iostream>
#include <memory>
#include <vector>

#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/environment/uniform_grid_environment.h"
#include "core/model_initializer.h"
#include "core/operation/mechanical_forces_op.h"
#include "core/operation/operation.h"
#include "core/param/param.h"
#include "core/real_t.h"
#include "core/simulation.h"

#include "cart_tumor.h"
#include "diffusion_thomas_algorithm.h"
#include "forces_tumor_cart.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "utils_aux.h"

namespace bdm {

int Simulate(int argc, const char** argv) {
  // Set simulation bounds
  auto set_param = [](Param* param) {
    param->random_seed = kSeed;  // Set a fixed random seed for reproducibility
    param->bound_space = Param::BoundSpaceMode::kTorus;  // Periodic boundary
    param->min_bound =
        -kBoundedSpaceLength /
        2.0;  // NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    param->max_bound =
        kBoundedSpaceLength /
        2.0;  // NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
              // Cube of 1000x1000x1000 centered at origin
    param->simulation_time_step = kDt;
  };

  Simulation simulation(argc, argv, set_param);
  auto* ctxt = simulation.GetExecutionContext();

  // Change Forces
  auto* scheduler = simulation.GetScheduler();

  auto* op = scheduler->GetOps("mechanical forces")[0];
  std::unique_ptr<InteractionVelocity> interaction_velocity =
      std::make_unique<InteractionVelocity>();
  op->GetImplementation<MechanicalForcesOp>()->SetInteractionForce(
      interaction_velocity.release());

  auto* env = dynamic_cast<UniformGridEnvironment*>(
      Simulation::GetActive()->GetEnvironment());
  // Fix the box length for the uniform grid environment
  env->SetBoxLength(gKLengthBoxMechanics);

  // Define Substances
  auto* rm = Simulation::GetActive()->GetResourceManager();

  // Oxygen
  // substance_id, name, diffusion_coefficient, decay_constant, resolution,
  // time_step
  std::unique_ptr<DiffusionThomasAlgorithm> oxygen_grid =
      std::make_unique<DiffusionThomasAlgorithm>(
          kOxygen, "oxygen",
          kDiffusionCoefficientOxygen,  // 100000 micrometers^2/minute
          kDecayConstantOxygen,         // 0.1 minutes^-1
          kResolutionGridSubstances, kDtSubstances,
          true);  // true indicates Dirichlet border conditions
  rm->AddContinuum(oxygen_grid.release());

  // Immunostimulatory Factor
  // substance_id, name, diffusion_coefficient, decay_constant, resolution
  std::unique_ptr<DiffusionThomasAlgorithm> immunostimulatory_factor_grid =
      std::make_unique<DiffusionThomasAlgorithm>(
          kImmunostimulatoryFactor, "immunostimulatory_factor",
          kDiffusionCoefficientImmunostimulatoryFactor,  // 1000
                                                         // micrometers^2/minute
          kDecayConstantImmunostimulatoryFactor,         // 0.016 minutes^-1
          kResolutionGridSubstances, kDtSubstances,
          false);  // false indicates Neumann border conditions
  rm->AddContinuum(immunostimulatory_factor_grid.release());

  // Boundary Conditions Dirichlet: simulating absorption or total loss at the
  // boundaries of the space.
  // Oxygen comming from the borders (capillary vessels)
  ModelInitializer::AddBoundaryConditions(
      kOxygen, BoundaryConditionType::kDirichlet,
      std::make_unique<ConstantBoundaryCondition>(
          kOxygenReferenceLevel));  // kOxygenReferenceLevel mmHg is the
                                    // physiological level of oxygen in tissues,
                                    // o2 saturation is 100% at this level

  // This is useless now but should be added this way in a future version of
  // BioDynaMo
  ModelInitializer::AddBoundaryConditions(
      kImmunostimulatoryFactor, BoundaryConditionType::kNeumann, nullptr);

  // Initialize oxygen voxels
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  ModelInitializer::InitializeSubstance(kOxygen, [](real_t /*x*/, real_t /*y*/,
                                                    real_t /*z*/) {
    return kInitialOxygenLevel;  // Set all voxels to kInitialOxygenLevel mmHg
  });

  // One spherical tumor of radius kInitialRadiusTumor in the center of the
  // simulation space
  const std::vector<Real3> positions =
      CreateSphereOfTumorCells(kInitialRadiusTumor);
  for (const auto& pos : positions) {
    std::unique_ptr<TumorCell> tumor_cell = std::make_unique<TumorCell>(pos);
    std::unique_ptr<StateControlGrowProliferate> state_control =
        std::make_unique<StateControlGrowProliferate>();
    tumor_cell->AddBehavior(state_control.release());
    ctxt->AddAgent(tumor_cell.release());
  }

  // OutputSummary operation
  std::unique_ptr<bdm::Operation> summary_op =
      std::make_unique<bdm::Operation>("OutputSummary", kOutputCsvInterval);
  std::unique_ptr<bdm::OutputSummary> output_summary =
      std::make_unique<bdm::OutputSummary>();
  summary_op->AddOperationImpl(bdm::kCpu, output_summary.release());
  scheduler->ScheduleOp(summary_op.release());

  // Run simulation
  // simulate kTotalMinutesToSimulate minutes including the last minute
  scheduler->Simulate(1 + kTotalMinutesToSimulate / kDt);
  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

int main(int argc, const char** argv) { return bdm::Simulate(argc, argv); }