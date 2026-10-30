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
#include "diffusion/cylinder_wall_boundary_condition.h"  // NOLINT(llvm-include-order)
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
#include "core/util/log.h"
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

namespace bdm {

namespace {


std::unique_ptr<Simulation> CreateSimulation(
    int argc, const char** argv, const std::function<void(Param*)>& set_param) {
  if (argv == nullptr) {
    return std::make_unique<Simulation>("CARTopiaX", set_param);
  }
  return std::make_unique<Simulation>(argc, argv, set_param);
}

}  // namespace

int RunSimulation(const SimParam& input_params, SimulationResult* result,
                  bool write_csv, bool verbose, bool allow_visualization,
                  int argc, const char** argv) {
  // Copy the parameters of this run
  std::unique_ptr<SimParam> custom_parameters =
      std::make_unique<SimParam>(input_params);

  // Keep a reference to the parameters before releasing
  const auto* sparam_ref = custom_parameters.get();

  // Transfer ownership to BioDynaMo parameter system
  Param::RegisterParamGroup(custom_parameters.release());

  // Set simulation parameters using lambda
  auto set_param = [sparam_ref, allow_visualization](Param* param) {
    // Set simulation bounds using the parameters
    param->random_seed = sparam_ref->seed;
    param->bound_space = sparam_ref->bound_space_toplogy;
    param->min_bound = -sparam_ref->bounded_space_length / kHalf;
    param->max_bound = sparam_ref->bounded_space_length / kHalf;
    param->simulation_time_step = sparam_ref->dt_step;
    param->statistics = sparam_ref->output_performance_statistics;
    // Disable the visualization export of bdm.toml and the command line
    // without editing the file (e.g. for calibration runs)
    param->export_visualization =
        param->export_visualization && allow_visualization;
  };

  auto simulation = CreateSimulation(argc, argv, set_param);
  const auto* sparam = simulation->GetParam()->Get<SimParam>();
  // Print parameters
  if (verbose) {
    sparam->PrintParams();
  }

  ExecutionContext* ctxt = simulation->GetExecutionContext();

  // Change Forces
  Scheduler* scheduler = simulation->GetScheduler();

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

  // Obtain the tumor shape from the parameters
  const TumorShape tumor_shape = sparam->tumor_shape;

  // Oxygen
  // substance_id, name, diffusion_coefficient, decay_constant, resolution,
  // time_step
  std::unique_ptr<DiffusionThomasAlgorithm> oxygen_grid =
      std::make_unique<DiffusionThomasAlgorithm>(
          kOxygen, "oxygen", sparam->diffusion_coefficient_oxygen,
          sparam->decay_constant_oxygen, sparam->resolution_grid_substances,
          sparam->dt_substances,
          /*dirichlet_border=*/true, sparam->diffuse_oxygen_on_z_axis);
  rm->AddContinuum(oxygen_grid.release());

  // Boundary Conditions Dirichlet: simulating absorption or total loss at the
  // boundaries of the space.
  // Oxygen comming from the borders (capillary vessels).
  // The boundary condition secretes oxygen wherever z falls within
  // [lateral_oxygen_production_min_z, lateral_oxygen_production_max_z],
  // regardless of which face of the domain is being evaluated (x, y, or z
  // faces). Depending on how these two parameters are set, this results in
  // either:
  //  - Oxygen produced only from a horizontal strip of the side walls, if
  //    lateral_oxygen_production_min_z, lateral_oxygen_production_max_z
  //    are set strictly inside the domain (floor and roof fall outside the
  //    range and produce no oxygen), or
  //  - Oxygen also produced from the floor and/or roof, if
  //    lateral_oxygen_production_min_z equals the domain's minimum z
  //    (floor) and/or lateral_oxygen_production_max_z equals the
  //    domain's maximum z (roof).
  ModelInitializer::AddBoundaryConditions(
      kOxygen, BoundaryConditionType::kDirichlet,
      // oxygen_reference_level mmHg is the physiological level of oxygen in
      // tissues, o2 saturation is 100% at this level
      std::make_unique<CylinderWallBoundaryCondition>(
          sparam->oxygen_reference_level,
          sparam->lateral_oxygen_production_min_z,
          sparam->lateral_oxygen_production_max_z));

  const real_t min_initial_z_substances = sparam->min_initial_z_substances;
  const real_t max_initial_z_substances = sparam->max_initial_z_substances;
  const real_t initial_oxygen_level = sparam->initial_oxygen_level;
  // Initialize oxygen voxels
  ModelInitializer::InitializeSubstance(
      kOxygen,
      [initial_oxygen_level, min_initial_z_substances,
       max_initial_z_substances](real_t /*x*/, real_t /*y*/, real_t z) {
        if (z <= max_initial_z_substances && z >= min_initial_z_substances) {
          return initial_oxygen_level;
        }
        return 0.0;
      });

  // Immunostimulatory Factor
  if (sparam->add_immunostimulatory_factor) {
    // substance_id, name, diffusion_coefficient, decay_constant, resolution
    std::unique_ptr<DiffusionThomasAlgorithm> immunostimulatory_factor_grid =
        std::make_unique<DiffusionThomasAlgorithm>(
            kImmunostimulatoryFactor, "immunostimulatory_factor",
            sparam->diffusion_coefficient_immunostimulatory_factor,
            sparam->decay_constant_immunostimulatory_factor,
            sparam->resolution_grid_substances, sparam->dt_substances,
            /*dirichlet_border=*/false,
            sparam->diffuse_immunostimulatory_factor_on_z_axis);
    rm->AddContinuum(immunostimulatory_factor_grid.release());

    // This is useless now but should be added this way in a future version of
    // BioDynaMo. Neumann boundary condition: the is not consumed nor produced
    // at the boundaries
    ModelInitializer::AddBoundaryConditions(
        kImmunostimulatoryFactor, BoundaryConditionType::kNeumann, nullptr);
  }
  // Glucose
  // substance_id, name, diffusion_coefficient, decay_constant, resolution,
  // time_step
  if (sparam->add_glucose) {
    std::unique_ptr<DiffusionThomasAlgorithm> glucose_grid =
        std::make_unique<DiffusionThomasAlgorithm>(
            kGlucose, "glucose", sparam->diffusion_coefficient_glucose,
            sparam->decay_constant_glucose, sparam->resolution_grid_substances,
            sparam->dt_substances,
            /*dirichlet_border=*/false, sparam->diffuse_glucose_on_z_axis);
    rm->AddContinuum(glucose_grid.release());

    // This is useless now but should be added this way in a future version of
    // BioDynaMo. Neumann boundary condition: the is not consumed nor produced
    // at the boundaries
    ModelInitializer::AddBoundaryConditions(
        kGlucose, BoundaryConditionType::kNeumann, nullptr);
    const real_t initial_glucose_level = sparam->initial_glucose_level;
    const real_t squared_radius_glucose_initialization =
        sparam->max_radius_glucose_initialization *
        sparam->max_radius_glucose_initialization;
    // Initialize oxygen voxels
    ModelInitializer::InitializeSubstance(
        kGlucose,
        [initial_glucose_level, min_initial_z_substances,
         max_initial_z_substances, tumor_shape,
         squared_radius_glucose_initialization](real_t x, real_t y, real_t z) {
          if (z <= max_initial_z_substances && z >= min_initial_z_substances) {
            if (tumor_shape == TumorShape::kSphere &&
                (x * x + y * y + z * z <=
                 squared_radius_glucose_initialization)) {
              // Its spherical and the voxel is within the maximun spherical
              // radius
              return initial_glucose_level;
            }
            if (tumor_shape == TumorShape::kCylinder &&
                (x * x + y * y <= squared_radius_glucose_initialization)) {
              // Its cylindrical and the voxel is within the maximun cylindrical
              // radius
              return initial_glucose_level;
            }
          }
          return 0.0;
        });
  }

  // Tumor cells initialization
  std::vector<Real3> positions;
  switch (tumor_shape) {
    case TumorShape::kSphere: {
      // One spherical tumor of radius initial_spherical_tumor_radius in the
      // center of the simulation space
      positions =
          CreateSphereOfTumorCells(sparam->initial_spherical_tumor_radius);
      break;
    }

    case TumorShape::kCylinder: {
      // One cylindrical tumor of radius initial_spherical_tumor_radius and
      // height cylindrical_tumor_height in the center of the simulation space
      positions = CreateCylinderOfTumorCells(
          sparam->cylindrical_tumor_radius, sparam->cylindrical_tumor_height,
          sparam->initial_number_of_cylindrical_tumor_cells);
      break;
    }

    default:
      Log::Error("RunSimulation",
                 "Unknown tumor shape, please use 'sphere' or 'cylinder'.");
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
  output_summary->SetOutputConfig(result, write_csv);
  summary_op->AddOperationImpl(bdm::kCpu, output_summary.release());
  scheduler->ScheduleOp(summary_op.release());

  // Run simulation
  if (verbose) {
    std::cout << "Running simulation..." << std::endl;
  }
  // simulate total_minutes_to_simulate minutes including the last minute
  scheduler->Simulate(1 +
                      static_cast<uint64_t>(sparam->total_minutes_to_simulate /
                                            sparam->dt_step));
  if (verbose) {
    std::cout << "Simulation completed successfully!" << std::endl;
  }
  const std::string output_dir = simulation->GetOutputDir();
  simulation.reset();
  std::error_code error;
  // Only removes the directory if it is empty
  std::filesystem::remove(output_dir, error);
  return 0;
}

int Simulate(int argc, const char** argv) {
  // Load parameters from JSON file or use default values
  SimParam custom_parameters;
  custom_parameters.LoadParams("params.json");
  return RunSimulation(custom_parameters, /*result=*/nullptr,
                       /*write_csv=*/true, /*verbose=*/true,
                       /*allow_visualization=*/true, argc, argv);
}

}  // namespace bdm

int main(int argc, const char** argv) { return bdm::Simulate(argc, argv); }
