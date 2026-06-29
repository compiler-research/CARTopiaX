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

#include "cart_tumor.h"
#include "agents/tumor_cell.h"
#include "diffusion/diffusion_thomas_algorithm.h"
#include "forces/forces_tumor_cart.h"
#include "params/hyperparams.h"
#include "utils/utils_aux.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/environment/uniform_grid_environment.h"
#include "core/model_initializer.h"
#include "core/operation/mechanical_forces_op.h"
#include "core/operation/operation.h"
#include "core/param/param.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include "core/scheduler.h"
#include "core/simulation.h"
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

namespace bdm {

int Simulate(int argc, const char** argv) {
  // Load parameters from JSON file or use default values
  std::unique_ptr<SimParam> custom_parameters = std::make_unique<SimParam>();
  custom_parameters->LoadParams("params.json");

  // Keep a reference to the parameters before releasing
  const auto* sparam_ref = custom_parameters.get();

  // Transfer ownership to BioDynaMo parameter system
  Param::RegisterParamGroup(custom_parameters.release());

  // Set simulation parameters using lambda
  auto set_param = [sparam_ref](Param* param) {
    // Set simulation bounds using the parameters
    param->random_seed = sparam_ref->seed;
    param->bound_space = Param::BoundSpaceMode::kTorus;
    param->min_bound = -sparam_ref->bounded_space_length / kHalf;
    param->max_bound = sparam_ref->bounded_space_length / kHalf;
    param->simulation_time_step = sparam_ref->dt_step;
    param->statistics = sparam_ref->output_performance_statistics;
  };

  Simulation simulation(argc, argv, set_param);
  const auto* sparam = simulation.GetParam()->Get<SimParam>();
  // Print parameters
  sparam->PrintParams();

  ExecutionContext* ctxt = simulation.GetExecutionContext();

  // Change Forces
  Scheduler* scheduler = simulation.GetScheduler();

  Operation* op = scheduler->GetOps("mechanical forces")[0];
  std::unique_ptr<InteractionVelocity> interaction_velocity =
      std::make_unique<InteractionVelocity>();
  op->GetImplementation<MechanicalForcesOp>()->SetInteractionForce(
      interaction_velocity.release());

  auto* env = dynamic_cast<UniformGridEnvironment*>(
      Simulation::GetActive()->GetEnvironment());
  // Fix the box length for the uniform grid environment
  env->SetBoxLength(sparam->length_box_mechanics);

  // Define Substances
  ResourceManager* rm = Simulation::GetActive()->GetResourceManager();

  // Oxygen
  // substance_id, name, diffusion_coefficient, decay_constant, resolution,
  // time_step
  std::unique_ptr<DiffusionThomasAlgorithm> oxygen_grid =
      std::make_unique<DiffusionThomasAlgorithm>(
          kOxygen, "oxygen", sparam->diffusion_coefficient_oxygen,
          sparam->decay_constant_oxygen, sparam->resolution_grid_substances,
          sparam->dt_substances,
          /*dirichlet_border=*/true);
  rm->AddContinuum(oxygen_grid.release());

  // Immunostimulatory Factor
  // substance_id, name, diffusion_coefficient, decay_constant, resolution
  std::unique_ptr<DiffusionThomasAlgorithm> immunostimulatory_factor_grid =
      std::make_unique<DiffusionThomasAlgorithm>(
          kImmunostimulatoryFactor, "immunostimulatory_factor",
          sparam->diffusion_coefficient_immunostimulatory_factor,
          sparam->decay_constant_immunostimulatory_factor,
          sparam->resolution_grid_substances, sparam->dt_substances,
          /*dirichlet_border=*/false);
  rm->AddContinuum(immunostimulatory_factor_grid.release());

  // Boundary Conditions Dirichlet: simulating absorption or total loss at the
  // boundaries of the space.
  // Oxygen comming from the borders (capillary vessels)
  ModelInitializer::AddBoundaryConditions(
      kOxygen, BoundaryConditionType::kDirichlet,
      // oxygen_reference_level mmHg is the physiological level of oxygen in
      // tissues, o2 saturation is 100% at this level
      std::make_unique<ConstantBoundaryCondition>(
          sparam->oxygen_reference_level));

  // This is useless now but should be added this way in a future version of
  // BioDynaMo
  ModelInitializer::AddBoundaryConditions(
      kImmunostimulatoryFactor, BoundaryConditionType::kNeumann, nullptr);

  // Initialize oxygen voxels
  ModelInitializer::InitializeSubstance(
      kOxygen, [sparam](real_t /*x*/, real_t /*y*/, real_t /*z*/) {
        // Set all voxels to initial_oxygen_level mmHg
        return sparam->initial_oxygen_level;
      });

  // Tumor cells initialization
  const std::string tumor_shape = sparam->tumor_shape;
  std::vector<Real3> positions;
  if (tumor_shape == "sphere") {
    // One spherical tumor of radius initial_spherical_tumor_radius in the center of the
    // simulation space
    positions = CreateSphereOfTumorCells(sparam->initial_spherical_tumor_radius);
  } else if (tumor_shape == "cylinder") {
    // One cylindrical tumor of radius initial_spherical_tumor_radius and height
    // cylindrical_tumor_height in the center of the simulation space
    positions = CreateCylinderOfTumorCells(
        sparam->cylindrical_tumor_radius, sparam->cylindrical_tumor_height,
        sparam->initial_number_of_cylindrical_tumor_cells);
  } else {
    Log::Error("Simulate", "Unknown tumor shape, please use 'sphere' or 'cylinder'.");
    // Exit with an error code
    return 1;
  }

  for (const auto& pos : positions) {
    std::unique_ptr<TumorCell> tumor_cell = std::make_unique<TumorCell>(pos);
    std::unique_ptr<StateControlGrowProliferate> state_control =
        std::make_unique<StateControlGrowProliferate>();
    tumor_cell->AddBehavior(state_control.release());
    ctxt->AddAgent(tumor_cell.release());
  }

  // Treatment administration operation
  std::unique_ptr<bdm::Operation> treatment_op =
      std::make_unique<bdm::Operation>("SpawnCart", sparam->steps_in_one_day);
  std::unique_ptr<bdm::SpawnCart> spawn_cart =
      std::make_unique<bdm::SpawnCart>();
  treatment_op->AddOperationImpl(bdm::kCpu, spawn_cart.release());
  scheduler->ScheduleOp(treatment_op.release());

  // OutputSummary operation
  std::unique_ptr<bdm::Operation> summary_op = std::make_unique<bdm::Operation>(
      "OutputSummary", sparam->output_csv_interval);
  std::unique_ptr<bdm::OutputSummary> output_summary =
      std::make_unique<bdm::OutputSummary>();
  summary_op->AddOperationImpl(bdm::kCpu, output_summary.release());
  scheduler->ScheduleOp(summary_op.release());

  // Run simulation
  std::cout << "Running simulation..." << std::endl;
  // simulate total_minutes_to_simulate minutes including the last minute
  scheduler->Simulate(1 +
                      static_cast<uint64_t>(sparam->total_minutes_to_simulate /
                                            sparam->dt_step));
  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

int main(int argc, const char** argv) { return bdm::Simulate(argc, argv); }
