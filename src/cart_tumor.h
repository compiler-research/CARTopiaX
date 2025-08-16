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



/**
 * @brief Enumeration of extracellular substances in the simulation
 * 
 * This enum defines the different types of substances that can diffuse
 * through the extracellular environment in the simulation.
 */
enum Substances { 
  /** @brief Immunostimulatory factor substance identifier */
  kImmunostimulatoryFactor, 
  /** @brief Oxygen substance identifier */
  kOxygen 
};

/**
 * @brief Main simulation function for CAR-T cell and tumor interaction
 * 
 * This function sets up and runs the complete simulation including:
 * - Simulation parameters and boundary conditions
 * - Mechanical forces and grid environment
 * - Diffusion grids for oxygen and immunostimulatory factors
 * - Initial tumor cell population
 * - Simulation execution and output
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return Exit code (0 for success)
 */
inline int Simulate(int argc, const char** argv) {



  /** @brief Configure simulation parameters
   * 
   * Lambda function to set up simulation parameters including:
   * - Random seed for reproducibility
   * - Boundary conditions (torus/periodic)
   * - Spatial bounds and simulation time step
   */
  auto set_param = [](Param* param) {
    param->random_seed = kSeed; /** @brief Set a fixed random seed for reproducibility */
    param->bound_space = Param::BoundSpaceMode::kTorus; /** @brief Periodic boundary conditions */
    param->min_bound = -kBoundedSpaceLength / 2;
    param->max_bound = kBoundedSpaceLength/2;  /** @brief Cube of 1000x1000x1000 centered at origin */
    param->simulation_time_step = kDt;
  };
  


  Simulation simulation(argc, argv, set_param);
  auto* ctxt = simulation.GetExecutionContext();

  /** @brief Configure mechanical forces
   * 
   * Change the default mechanical forces to use custom interaction velocity forces
   * and set up the uniform grid environment with specified box length.
   */
  auto* scheduler = simulation.GetScheduler();

  auto* op = scheduler->GetOps("mechanical forces")[0];
  op->GetImplementation<MechanicalForcesOp>()->SetInteractionForce(new InteractionVelocity());

  auto* env = dynamic_cast<UniformGridEnvironment*>(Simulation::GetActive()->GetEnvironment());
  env->SetBoxLength(kLengthBoxMechanics); /** @brief Fix the box length for the uniform grid environment */

  /** @name Substance Definition and Configuration
   *  @brief Setup of diffusion grids for extracellular substances
   *  @{
   */
  auto* rm = Simulation::GetActive()->GetResourceManager();

  /** @brief Oxygen diffusion grid setup
   * 
   * Creates a diffusion grid for oxygen with Thomas algorithm solver.
   * Parameters: substance_id, name, diffusion_coefficient, decay_constant, resolution, time_step
   * Uses Dirichlet boundary conditions to simulate oxygen supply from capillary vessels.
   */
  auto* oxygen_grid = new DiffusionThomasAlgorithm(
      kOxygen, "oxygen",
      kDiffusionCoefficientOxygen,/** @brief 100000 micrometers^2/minute */
      kDecayConstantOxygen, /** @brief 0.1 minutes^-1 */
      kResolutionGridSubstances,
      kDtSubstances,
      true); /** @brief true indicates Dirichlet border conditions */
  rm->AddContinuum(oxygen_grid);

  /** @brief Immunostimulatory factor diffusion grid setup
   * 
   * Creates a diffusion grid for immunostimulatory factors with Thomas algorithm solver.
   * Parameters: substance_id, name, diffusion_coefficient, decay_constant, resolution
   * Uses Neumann boundary conditions (no flux across boundaries).
   */
  auto* immunostimulatory_factor_grid = new DiffusionThomasAlgorithm(
      kImmunostimulatoryFactor, "immunostimulatory_factor",
      kDiffusionCoefficientImmunostimulatoryFactor, /** @brief 1000 micrometers^2/minute */
      kDecayConstantImmunostimulatoryFactor, /** @brief 0.016 minutes^-1 */
      kResolutionGridSubstances,
      kDtSubstances,
    false); /** @brief false indicates Neumann border conditions */
  rm->AddContinuum(immunostimulatory_factor_grid);


  /** @brief Boundary conditions setup
   * 
   * Dirichlet boundary conditions simulate absorption or total loss at the boundaries.
   * Oxygen comes from the borders (simulating capillary vessels).
   */
  ModelInitializer::AddBoundaryConditions(
    kOxygen, BoundaryConditionType::kDirichlet,
    std::make_unique<ConstantBoundaryCondition>(kOxygenReferenceLevel));/** @brief kOxygenReferenceLevel mmHg is the physiological level of oxygen in tissues, O2 saturation is 100% at this level */

  /** @brief Neumann boundary conditions for immunostimulatory factor
   * 
   * This is currently not used but should be added this way in a future version of BioDynaMo
   */
  ModelInitializer::AddBoundaryConditions(
      kImmunostimulatoryFactor, BoundaryConditionType::kNeumann, nullptr);

  /** @brief Initialize oxygen concentration throughout the simulation space */
  ModelInitializer::InitializeSubstance(kOxygen, [](real_t x, real_t y, real_t z) {
    return kInitialOxygenLevel; /** @brief Set all voxels to kInitialOxygenLevel mmHg */
  });

  /** @} */ // end of Substance Definition and Configuration group

  /** @name Initial Cell Population Setup
   *  @brief Creation of initial tumor cell population
   *  @{
   */
  
  /** @brief Create spherical tumor of specified radius
   * 
   * Generate positions for tumor cells arranged in a sphere of radius kInitialRadiusTumor
   * centered at the origin of the simulation space.
   */
  std::vector<Real3> positions=CreateSphereOfTumorCells(kInitialRadiusTumor);/** @brief positions of the cancer cells */

  /** @brief Create and initialize tumor cells */
  for (const auto& pos : positions) {
    TumorCell* tumor_cell = new TumorCell(pos);
    tumor_cell->AddBehavior(new StateControlGrowProliferate());
    ctxt->AddAgent(tumor_cell);
  }

  /** @brief Debug CAR-T cell creation (commented out)
   * 
   * Uncomment to add a single CAR-T cell at origin for debugging purposes
   */
  // //debug
  // CartCell* cart_cell = new CartCell({0.,0.,0.});
  // cart_cell->AddBehavior(new StateControlCart());
  // ctxt->AddAgent(cart_cell);

  /** @} */ // end of Initial Cell Population Setup group




  /** @name Output Operations Setup
   *  @brief Configuration of simulation output and data collection
   *  @{
   */
  
  /** @brief Setup output summary operation
   * 
   * Creates and schedules an operation to output CSV files at regular intervals
   * for data analysis and visualization.
   */
  auto* summary_op = new bdm::Operation("OutputSummary");
  summary_op->frequency_ = kOutputCsvInterval; /** @brief Set the interval for outputting CSV files */
  summary_op->AddOperationImpl(bdm::kCpu, new bdm::OutputSummary());
  scheduler->ScheduleOp(summary_op);

  /** @} */ // end of Output Operations Setup group



  /** @name Simulation Execution
   *  @brief Run the main simulation loop
   *  @{
   */
  
  /** @brief Execute the simulation
   * 
   * Run the simulation for the specified duration. The simulation runs for
   * kTotalMinutesToSimulate minutes including the final minute (hence 1+).
   */
  scheduler->Simulate(1+kTotalMinutesToSimulate/kDt);/** @brief simulate kTotalMinutesToSimulate minutes including the last minute */
  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;

  /** @} */ // end of Simulation Execution group
}

}  // namespace bdm

#endif  // CART_TUMOR_H_
