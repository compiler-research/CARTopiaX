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

#ifndef TUMOR_CELL_H_
#define TUMOR_CELL_H_

#include "biodynamo.h"
#include "core/util/log.h"
#include "core/util/root.h"
#include "hyperparams.h"
#include "utils_aux.h"

namespace bdm {

/**
 * @brief Enumeration representing the different states of a tumor cell
 * 
 * This enum class defines the various states a tumor cell can be in during its lifecycle,
 * including different proliferation phases based on Ki67 protein expression and various
 * death pathways (necrosis and apoptosis).
 */
enum class TumorCellState : int {
  /** @brief Living cell state - cell is alive and can potentially proliferate */
  kAlive=0,

  /** @brief Necrotic swelling phase
   *  
   *  The cell loses membrane integrity and starts absorbing fluid, swelling abnormally
   *  in volume before rupture. This is the first phase of necrotic cell death.
   */
  kNecroticSwelling = 1,
  
  /** @brief Necrotic lysed phase
   *  
   *  The cell membrane breaks apart, releasing its contents. The cell is now considered
   *  dead and will be removed from the simulation after a defined time.
   */
  kNecroticLysed = 2,
    
  /** @brief Apoptotic phase
   *  
   *  The cell is undergoing programmed cell death characterized by cell shrinkage.
   *  This is a controlled form of cell death.
   */
  kApoptotic=3
};

/**
 * @brief Tumor cell class implementation
 * 
 * This class represents a tumor cell in the simulation with capabilities for:
 * - Different cellular states (alive, necrotic, apoptotic)
 * - Volume dynamics with exponential relaxation
 * - Oxygen consumption and immunostimulatory factor secretion
 * - Oncoprotein expression levels
 * - Interactions with CAR-T cells
 */
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);

  public:
  /** @brief Default constructor */
  TumorCell() {}
  
  /** @brief Constructor with position parameter
   *  @param position Initial 3D position of the cell
   */
  explicit TumorCell(const Real3& position);
  
  /** @brief Virtual destructor */
  virtual ~TumorCell() {}

  /** @brief Called when a new agent is created (e.g., after cell division)
   *  @param event The new agent event containing initialization data
   */
  void Initialize(const NewAgentEvent& event) override;

  /** @name State Management
   *  @brief Methods for managing tumor cell state
   *  @{
   */

  /** @brief Set the current state of the tumor cell
   *  @param state The new state to set
   */
  void SetState(TumorCellState state) { state_ = state; }
  
  /** @brief Get the current state of the tumor cell
   *  @return The current cell state
   */
  TumorCellState GetState() const { return state_; }

  /** @brief Set the timer for tracking time in current state
   *  @param timer_state Timer value in minutes
   */
  void SetTimerState(int timer_state) { timer_state_ = timer_state; }
  
  /** @brief Get the timer for tracking time in current state
   *  @return Timer value in minutes
   */
  int GetTimerState() const { return timer_state_; }

  /** @} */ // end of State Management group

  /** @name Oncoprotein Management
   *  @brief Methods for managing oncoprotein expression levels
   *  @{
   */

  /** @brief Set the oncoprotein expression level
   *  @param level The oncoprotein level to set
   */
  void SetOncoproteineLevel(real_t level);
  
  /** @brief Get the oncoprotein expression level
   *  @return The current oncoprotein level
   */
  real_t GetOncoproteineLevel() const { return oncoproteine_level_; }

  /** @} */ // end of Oncoprotein Management group  

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

  /** @name Transformation and Cell Type
   *  @brief Methods for managing transformation rates and cell types
   *  @{
   */

  /** @brief Set the transformation random rate for state transitions
   * 
   *  This rate affects the probability of transitioning between states and
   *  depends on the individual cell. It remains constant during the cell's lifetime.
   */
  void SetTransformationRandomRate();
  
  /** @brief Get the transformation random rate
   *  @return The current transformation random rate
   */
  real_t GetTransformationRandomRate() const { return transformation_random_rate_; }

  /** @brief Set the cell type based on oncoprotein level
   *  @param type Cell type (1-4: 1 is most mutated/proliferative, 4 is least aggressive; 5 means dead)
   */
  void SetType(int type) { type_ = type; }
  
  /** @brief Get the cell type
   *  @return The current cell type
   */
  int GetType() const { return type_; }

  /** @name CAR-T Cell Interaction
   *  @brief Methods for managing attachment to CAR-T cells
   *  @{
   */

  /** @brief Set whether the cell is attached to a CAR-T cell
   *  @param attached True if attached, false otherwise
   */
  void SetAttachedToCart(bool attached) { attached_to_cart_ = attached; }
  
  /** @brief Check if the cell is attached to a CAR-T cell
   *  @return True if attached to a CAR-T cell, false otherwise
   */
  bool IsAttachedToCart() const { return attached_to_cart_; }

  /** @} */ // end of CAR-T Cell Interaction group

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

  /** @} */ // end of Movement and Velocity group

  /** @name Biochemical Properties
   *  @brief Methods for managing substance consumption and secretion
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

  /** @brief Get the immunostimulatory factor secretion rate
   *  @return The current immunostimulatory factor secretion rate
   */
  real_t GetImmunostimulatoryFactorSecretionRate() const { return immunostimulatory_factor_secretion_rate_; }
  
  /** @brief Set the immunostimulatory factor secretion rate
   *  @param rate The immunostimulatory factor secretion rate to set
   */
  void SetImmunostimulatoryFactorSecretionRate(real_t rate) { immunostimulatory_factor_secretion_rate_ = rate; }

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
  *  @brief Private attributes of the tumor cell
  *  @{
  */
 private:
  /** @brief Current state of the tumor cell */
  TumorCellState state_;
  
  /** @brief Timer to track time in the current state (in minutes) */
  int timer_state_;
  
  /** @brief Pointer to the oxygen diffusion grid */
  DiffusionGrid* oxygen_dgrid_;
  
  /** @brief Pointer to the immunostimulatory factor diffusion grid */
  DiffusionGrid* immunostimulatory_factor_dgrid_;
  
  /** @brief Level of oncoprotein expression */
  real_t oncoproteine_level_;
  
  /** @brief Transition random rate between states
   * 
   *  Affects the probability of transitioning and depends on the individual cell.
   *  This rate is kept constant during the cell's lifetime.
   */
  real_t transformation_random_rate_;
  
  /** @brief Flag indicating if the cell is attached to a CAR-T cell */
  bool attached_to_cart_;
  
  /** @brief Fluid fraction of the cell volume */
  real_t fluid_fraction_;
  
  /** @brief Volume of the nucleus */
  real_t nuclear_volume_;
  
  /** @brief Target cytoplasm solid volume for exponential relaxation
   *  
   *  Used for growing or shrinking tumor cells. The volume change follows
   *  an exponential relaxation equation toward this target volume.
   */
  real_t target_cytoplasm_solid_;
  
  /** @brief Target nucleus solid volume for exponential relaxation */
  real_t target_nucleus_solid_;
  
  /** @brief Target fluid fraction for exponential relaxation */
  real_t target_fraction_fluid_;
  
  /** @brief Target relation between cytoplasm and nucleus volumes */
  real_t target_relation_cytoplasm_nucleus_;
  
  /** @brief Cell type according to oncoprotein level
   * 
   *  Types 1-4: 1 is the most mutated and proliferative type, 4 is the least aggressive.
   *  Type 5 means the cell is dead.
   */
  int type_;
  
  /** @brief Velocity of the cell in the previous time step */
  Real3 older_velocity_;
  
  /** @brief Rate of oxygen consumption by the cell */
  real_t oxygen_consumption_rate_;
  
  /** @brief Rate of immunostimulatory factor secretion by the cell */
  real_t immunostimulatory_factor_secretion_rate_;
  
  /** @name Consumption/Secretion Constants
   *  @brief Constants for consumption/secretion differential equation solutions
   *  @{
   */
  
  /** @brief Constant 1 for oxygen consumption/secretion differential equation solution */
  real_t constant1_oxygen_;
  
  /** @brief Constant 2 for oxygen consumption/secretion differential equation solution */
  real_t constant2_oxygen_;
  
  /** @brief Constant 1 for immunostimulatory factor consumption/secretion differential equation solution */
  real_t constant1_immunostimulatory_factor_;
  
  /** @brief Constant 2 for immunostimulatory factor consumption/secretion differential equation solution */
  real_t constant2_immunostimulatory_factor_;
  
  /** @} */ // end of Consumption/Secretion Constants group

  /** @} */ // end of Private Member Variables group
};

/**
 * @brief Behavior class for controlling tumor cell state transitions and growth
 * 
 * This behavior handles the state control logic for tumor cells, managing
 * transitions between different cell states, growth, proliferation, and death
 * processes. It includes logic for determining when cells should enter necrosis
 * based on oxygen levels and other environmental factors.
 */
struct StateControlGrowProliferate : public Behavior {
  BDM_BEHAVIOR_HEADER(StateControlGrowProliferate, Behavior, 1);

  /** @brief Default constructor
   *  Calls AlwaysCopyToNew() to ensure the behavior is copied to new cells
   */
  StateControlGrowProliferate() { AlwaysCopyToNew(); }
  
  /** @brief Virtual destructor */
  virtual ~StateControlGrowProliferate() {}

  /** @brief Execute the state control and growth behavior
   *  @param agent Pointer to the agent (cell) on which to apply the behavior
   */
  void Run(Agent* agent) override;

  private:
  /** @brief Compute the probability of the cell entering necrosis
   * 
   *  Determines whether a cell should enter necrosis based on oxygen levels
   *  and other cellular conditions.
   * 
   *  @param oxygen_level Current oxygen concentration at the cell's location
   *  @param cell Pointer to the tumor cell being evaluated
   *  @return True if the cell should enter necrosis, false otherwise
   */
  bool ShouldEnterNecrosis(real_t oxygen_level,TumorCell* cell) const;
};

}  // namespace bdm

#endif  // TUMOR_CELL_H_