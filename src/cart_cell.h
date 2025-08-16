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

#ifndef CART_CELL_H_
#define CART_CELL_H_

#include "biodynamo.h"
#include "core/util/log.h"
#include "core/util/root.h"
#include "utils_aux.h"
#include "hyperparams.h"
#include "tumor_cell.h"

namespace bdm {


/**
 * @brief Enumeration defining the possible states of a CAR-T cell
 * 
 * This enum class represents the different states that a CAR-T cell can be in
 * during its lifecycle in the simulation.
 */
enum class CartCellState : int {
  /** @brief Living cell state - the cell is alive and functioning normally */
  kAlive=0,

  /** @brief Apoptotic phase - the cell is undergoing programmed cell death
   *  characterized by cell shrinkage and controlled death
   */
  kApoptotic=1
};

/**
 * @brief CAR-T cell class implementation
 * 
 * This class represents a CAR-T (Chimeric Antigen Receptor T-cell) in the simulation.
 * It inherits from the base Cell class and includes specific behaviors and properties
 * related to CAR-T cell biology, including states, volume dynamics, and interactions
 * with tumor cells.
 */
class CartCell : public Cell {
  BDM_AGENT_HEADER(CartCell, Cell, 1);

  public:
  /** @brief Default constructor */
  CartCell() {}
  
  /** @brief Constructor with position parameter
   *  @param position Initial 3D position of the cell
   */
  explicit CartCell(const Real3& position);
  
  /** @brief Virtual destructor */
  virtual ~CartCell() {}

  /** @name State Management
   *  @brief Methods for managing cell state
   *  @{
   */
  
  /** @brief Set the current state of the CAR-T cell
   *  @param state The new state to set
   */
  void SetState(CartCellState state) { state_ = state; }
  
  /** @brief Get the current state of the CAR-T cell
   *  @return The current cell state
   */
  CartCellState GetState() const { return state_; }

  /** @brief Set the timer for tracking time in current state
   *  @param timer_state Timer value in minutes
   */
  void SetTimerState(int timer_state) { timer_state_ = timer_state; }
  
  /** @brief Get the timer for tracking time in current state
   *  @return Timer value in minutes
   */
  int GetTimerState() const { return timer_state_; }
  
  /** @} */ // end of State Management group

  /** @name Volume and Physical Properties
   *  @brief Methods for managing cell volume and physical characteristics
   *  @{
   */

  /** @brief Set the fluid fraction of the cell
   *  @param fluid_fraction The fluid fraction value
   */
  void SetFluidFraction(real_t fluid_fraction) { fluid_fraction_ = fluid_fraction; }
  
  /** @brief Get the fluid fraction of the cell
   *  @return The current fluid fraction
   */
  real_t GetFluidFraction() const { return fluid_fraction_; }

  /** @brief Set the nuclear volume
   *  @param nuclear_volume The nuclear volume value
   */
  void SetNuclearVolume(real_t nuclear_volume) { nuclear_volume_ = nuclear_volume; }
  
  /** @brief Get the nuclear volume
   *  @return The current nuclear volume
   */
  real_t GetNuclearVolume() const { return nuclear_volume_; }

  /** @brief Set the target cytoplasm solid volume
   *  @param target_cytoplasm_solid The target cytoplasm solid volume
   */
  void SetTargetCytoplasmSolid(real_t target_cytoplasm_solid) { target_cytoplasm_solid_ = target_cytoplasm_solid; }
  
  /** @brief Get the target cytoplasm solid volume
   *  @return The target cytoplasm solid volume
   */
  real_t GetTargetCytoplasmSolid() const { return target_cytoplasm_solid_; }

  /** @brief Set the target nucleus solid volume
   *  @param target_nucleus_solid The target nucleus solid volume
   */
  void SetTargetNucleusSolid(real_t target_nucleus_solid) { target_nucleus_solid_ = target_nucleus_solid; }
  
  /** @brief Get the target nucleus solid volume
   *  @return The target nucleus solid volume
   */
  real_t GetTargetNucleusSolid() const { return target_nucleus_solid_; }  

  /** @brief Set the target fraction of fluid
   *  @param target_fraction_fluid The target fluid fraction
   */
  void SetTargetFractionFluid(real_t target_fraction_fluid) { target_fraction_fluid_ = target_fraction_fluid; }
  
  /** @brief Get the target fraction of fluid
   *  @return The target fluid fraction
   */
  real_t GetTargetFractionFluid() const { return target_fraction_fluid_; }  

  /** @brief Set the target relation between cytoplasm and nucleus
   *  @param target_relation_cytoplasm_nucleus The target relation value
   */
  void SetTargetRelationCytoplasmNucleus(real_t target_relation_cytoplasm_nucleus) { target_relation_cytoplasm_nucleus_ = target_relation_cytoplasm_nucleus; }
  
  /** @brief Get the target relation between cytoplasm and nucleus
   *  @return The target relation value
   */
  real_t GetTargetRelationCytoplasmNucleus() const { return target_relation_cytoplasm_nucleus_; }

  /** @} */ // end of Volume and Physical Properties group

  /** @name Tumor Cell Attachment
   *  @brief Methods for managing attachment to tumor cells
   *  @{
   */

  /** @brief Set whether the cell is attached to a tumor cell
   *  @param attached True if attached, false otherwise
   */
  void SetAttachedToTumorCell(bool attached) { attached_to_tumor_cell_ = attached; }
  
  /** @brief Check if the cell is attached to a tumor cell
   *  @return True if attached to a tumor cell, false otherwise
   */
  bool IsAttachedToTumorCell() const { return attached_to_tumor_cell_; }

  /** @brief Get the attached tumor cell
   *  @return Pointer to the attached tumor cell, or nullptr if not attached
   */
  TumorCell* GetAttachedCell() const { return attached_cell_; }
  
  /** @brief Set the attached tumor cell
   *  @param cell Pointer to the tumor cell to attach
   */
  void SetAttachedCell(TumorCell* cell) { attached_cell_ = cell; }

  /** @} */ // end of Tumor Cell Attachment group

  /** @name Movement and Velocity
   *  @brief Methods for managing cell movement and velocity
   *  @{
   */

  /** @brief Get the velocity from the previous time step
   *  @return The velocity vector from the previous step
   */
  Real3 GetOlderVelocity() const { return older_velocity_; }
  
  /** @brief Set the velocity from the previous time step
   *  @param velocity The velocity vector to set
   */
  void SetOlderVelocity(const Real3& velocity) { older_velocity_ = velocity; }

  /** @brief Check whether the cell moves by its own
   *  @return True if the cell can move independently, false otherwise
   */
  bool DoesCellMove();

  /** @} */ // end of Movement and Velocity group

  /** @name Biochemical Properties
   *  @brief Methods for managing oxygen consumption and cell lifetime
   *  @{
   */

  /** @brief Get the oxygen consumption rate
   *  @return The current oxygen consumption rate
   */
  real_t GetOxygenConsumptionRate() const { return oxygen_consumption_rate_; }
  
  /** @brief Set the oxygen consumption rate
   *  @param rate The oxygen consumption rate to set
   */
  void SetOxygenConsumptionRate(real_t rate) { oxygen_consumption_rate_ = rate; }

  /** @brief Get the current live time
   *  @return The current time until apoptosis
   */
  real_t GetCurrentLiveTime() const { return current_live_time_; }
  
  /** @brief Set the current live time
   *  @param time The current live time to set
   */
  void SetCurrentLiveTime(real_t time) { current_live_time_ = time; }

  /** @} */ // end of Biochemical Properties group

  /** @name Volume Calculations
   *  @brief Methods for volume calculations
   *  @{
   */

  /** @brief Calculate the target total volume of the cell
   *  @return The target total volume
   */
  real_t GetTargetTotalVolume();

  /** @} */ // end of Volume Calculations group

  /** @name Diffusion Grids
   *  @brief Methods for accessing diffusion grids
   *  @{
   */

  /** @brief Get the diffusion grid for oxygen
   *  @return Pointer to the oxygen diffusion grid
   */
  DiffusionGrid* GetOxygenDiffusionGrid() const { return oxygen_dgrid_; }
  
  /** @brief Get the diffusion grid for immunostimulatory factors
   *  @return Pointer to the immunostimulatory factor diffusion grid
   */
  DiffusionGrid* GetImmunostimulatoryFactorDiffusionGrid() const { return immunostimulatory_factor_dgrid_; }

  /** @} */ // end of Diffusion Grids group

  /** @name Core Simulation Methods
   *  @brief Core methods for cell simulation and behavior
   *  @{
   */

  /** @brief Change volume using exponential relaxation equation
   * 
   * This method explicitly solves the system of exponential relaxation differential
   * equations using a discrete update step. It is used to grow or shrink the volume
   * (and proportions) smoothly toward a desired target volume over time. The relaxation
   * rate controls the speed of convergence and dt=1 (the time_step).
   * 
   * @param relaxation_rate_cytoplasm Relaxation rate for cytoplasm volume changes
   * @param relaxation_rate_nucleus Relaxation rate for nucleus volume changes
   * @param relaxation_rate_fluid Relaxation rate for fluid volume changes
   */
  void ChangeVolumeExponentialRelaxationEquation(real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus, real_t relaxation_rate_fluid);

  /** @brief Calculate displacement of the cell
   * 
   * Computes the displacement of the cell based on interaction forces.
   * 
   * @param force Pointer to the interaction force object
   * @param squared_radius The squared radius of the cell
   * @param dt The time step for the simulation
   * @return The calculated displacement vector
   */
  Real3 CalculateDisplacement(const InteractionForce* force,
                              real_t squared_radius, real_t dt) override;

  /** @brief Consume or secrete substances
   * 
   * Computes new oxygen or immunostimulatory factor concentration after
   * consumption or secretion by the cell.
   * 
   * @param substance_id The ID of the substance (oxygen or immunostimulatory factor)
   * @param old_concentration The previous concentration of the substance
   * @return The new concentration after consumption/secretion
   */
  real_t ConsumeSecreteSubstance(int substance_id, real_t old_concentration);

  /** @brief Compute constants for consumption and secretion
   * 
   * Updates constants after the cell's change of volume or quantities.
   * These constants are used in the consumption/secretion differential equations.
   */
  void ComputeConstantsConsumptionSecretion();

  /** @} */ // end of Core Simulation Methods group

 /** @name Private Member Variables
  *  @brief Private attributes of the CAR-T cell
  *  @{
  */
 private:
  /** @brief Current state of the CAR-T cell */
  CartCellState state_;
  
  /** @brief Timer to track time in the current state (in minutes)
   *  Used for apoptotic state timing
   */
  int timer_state_;
  
  /** @brief Pointer to the oxygen diffusion grid */
  DiffusionGrid* oxygen_dgrid_;
  
  /** @brief Pointer to the immunostimulatory factor diffusion grid */
  DiffusionGrid* immunostimulatory_factor_dgrid_;
  
  /** @brief Flag indicating if the cell is attached to a tumor cell */
  bool attached_to_tumor_cell_;
  
  /** @brief Current time until apoptosis */
  real_t current_live_time_;
  
  /** @brief Fluid fraction of the cell volume */
  real_t fluid_fraction_;
  
  /** @brief Volume of the nucleus */
  real_t nuclear_volume_;
  
  /** @brief Target cytoplasm solid volume for exponential relaxation
   *  Used during volume changes following exponential relaxation equation
   */
  real_t target_cytoplasm_solid_;
  
  /** @brief Target nucleus solid volume for exponential relaxation */
  real_t target_nucleus_solid_;
  
  /** @brief Target fluid fraction for exponential relaxation */
  real_t target_fraction_fluid_;
  
  /** @brief Target relation between cytoplasm and nucleus volumes */
  real_t target_relation_cytoplasm_nucleus_;
  
  /** @brief Velocity of the cell in the previous time step */
  Real3 older_velocity_;
  
  /** @brief Rate of oxygen consumption by the cell */
  real_t oxygen_consumption_rate_;
  
  /** @brief Rate of immunostimulatory factor secretion by the cell */
  real_t immunostimulatory_factor_secretion_rate_;
  
  /** @brief Constant 1 for oxygen consumption/secretion differential equation solution */
  real_t constant1_oxygen_;
  
  /** @brief Constant 2 for oxygen consumption/secretion differential equation solution */
  real_t constant2_oxygen_;
  
  /** @brief Pointer to the attached tumor cell */
  TumorCell* attached_cell_;

  /** @} */ // end of Private Member Variables group
};

/**
 * @brief Behavior class for controlling CAR-T cell state transitions
 * 
 * This behavior handles the state control logic for CAR-T cells, managing
 * transitions between different cell states such as alive and apoptotic phases.
 * It inherits from the base Behavior class and implements the Run method to
 * execute the state control logic during simulation steps.
 */
struct StateControlCart : public Behavior {
  BDM_BEHAVIOR_HEADER(StateControlCart, Behavior, 1);

  /** @brief Default constructor
   *  Calls AlwaysCopyToNew() to ensure the behavior is copied to new cells
   */
  StateControlCart() { AlwaysCopyToNew(); }
  
  /** @brief Virtual destructor */
  virtual ~StateControlCart() {}

  /** @brief Execute the state control behavior
   *  @param agent Pointer to the agent (cell) on which to apply the behavior
   */
  void Run(Agent* agent) override;
};

}  // namespace bdm

#endif  // CART_CELL_H_
