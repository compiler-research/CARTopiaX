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

#ifndef TUMOR_CELL_H_
#define TUMOR_CELL_H_

#include "core/util/root.h"
#include "core/agent/agent.h"
#include "core/agent/agent.h"  // Ensure bdm::Agent is directly included
#include "core/agent/cell.h"
#include "core/behavior/behavior.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/agent/new_agent_event.h"
#include "core/real_t.h"
#include "core/interaction_force.h"
#include "core/resource_manager.h"
#include "core/scheduler.h"
#include "core/container/math_array.h"

namespace bdm {

/// Enumeration representing the different states of a tumor cell
///
/// This enum class defines the various states a tumor cell can be in during its
/// lifecycle, and various death pathways (necrosis and apoptosis).
enum class TumorCellState : int {
  kAlive =
      0,  ///< Living cell state - cell is alive and can potentially proliferate

  kNecroticSwelling =
      1,  ///< Necrotic swelling phase: The cell loses membrane integrity and
          ///< starts absorbing fluid, swelling abnormally, in volume before
          ///< rupture. This is the first phase of necrotic cell death.
  kNecroticLysed = 2,  ///< Necrotic lysed phase: The cell membrane breaks
                       ///< apart, releasing its contents. The cell will be
                       ///< removed from the simulation after a defined time.

  kApoptotic = 3  ///< Apoptotic phase: The cell is undergoing programmed cell
                  ///< death characterized by cell shrinkage. This is a
                  ///< controlled form of cell death.
};

/// Tumor cell class implementation
///
/// This class represents a cancer cell that forms a heterogeneous tumor in the
/// simulation. The class includes capabilities for:
/// - Different cellular states (alive, necrotic, apoptotic)
/// - Volume dynamics with exponential relaxation
/// - Cell division for tumor proliferation
/// - Oxygen consumption and immunostimulatory factor secretion
/// - Displacement computation applying pushing/adhesive forces between cells
/// - Oncoprotein expression levels
/// - Interactions with CAR-T cells

class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);

 public:
  TumorCell()
      : state_(TumorCellState::kAlive),
        timer_state_(0),
        oxygen_dgrid_(nullptr),
        immunostimulatory_factor_dgrid_(nullptr),
        oncoproteine_level_(0.0),
        transformation_random_rate_(0.0),
        attached_to_cart_(false),
        fluid_fraction_(0.0),
        nuclear_volume_(0.0),
        target_cytoplasm_solid_(0.0),
        target_nucleus_solid_(0.0),
        target_fraction_fluid_(0.0),
        target_relation_cytoplasm_nucleus_(0.0),
        type_(0),
        oxygen_consumption_rate_(0.0),
        immunostimulatory_factor_secretion_rate_(0.0),
        constant1_oxygen_(0.0),
        constant2_oxygen_(0.0),
        constant1_immunostimulatory_factor_(0.0),
        constant2_immunostimulatory_factor_(0.0) {}

  explicit TumorCell(const Real3& position);

  // Special member functions
  TumorCell(const TumorCell&) = default;
  TumorCell(TumorCell&&) = default;

  ~TumorCell() override = default;

  /// Called when a new agent is created (after cell division)
  /// @param event The new agent event containing initialization data
  void Initialize(const NewAgentEvent& event) override;

  /// Getters and Setters
  void SetState(TumorCellState state) { state_ = state; }
  TumorCellState GetState() const { return state_; }

  void SetTimerState(int timer_state) { timer_state_ = timer_state; }
  int GetTimerState() const { return timer_state_; }

  void SetOncoproteineLevel(real_t level);
  real_t GetOncoproteineLevel() const { return oncoproteine_level_; }

  void SetFluidFraction(real_t fluid_fraction) {
    fluid_fraction_ = fluid_fraction;
  }
  real_t GetFluidFraction() const { return fluid_fraction_; }

  void SetNuclearVolume(real_t nuclear_volume) {
    nuclear_volume_ = nuclear_volume;
  }
  real_t GetNuclearVolume() const { return nuclear_volume_; }

  void SetTargetCytoplasmSolid(real_t target_cytoplasm_solid) {
    target_cytoplasm_solid_ = target_cytoplasm_solid;
  }
  real_t GetTargetCytoplasmSolid() const { return target_cytoplasm_solid_; }

  void SetTargetNucleusSolid(real_t target_nucleus_solid) {
    target_nucleus_solid_ = target_nucleus_solid;
  }
  real_t GetTargetNucleusSolid() const { return target_nucleus_solid_; }

  void SetTargetFractionFluid(real_t target_fraction_fluid) {
    target_fraction_fluid_ = target_fraction_fluid;
  }
  real_t GetTargetFractionFluid() const { return target_fraction_fluid_; }

  void SetTargetRelationCytoplasmNucleus(
      real_t target_relation_cytoplasm_nucleus) {
    target_relation_cytoplasm_nucleus_ = target_relation_cytoplasm_nucleus;
  }
  real_t GetTargetRelationCytoplasmNucleus() const {
    return target_relation_cytoplasm_nucleus_;
  }

  void SetTransformationRandomRate();
  real_t GetTransformationRandomRate() const {
    return transformation_random_rate_;
  }

  void SetAttachedToCart(bool attached) { attached_to_cart_ = attached; }
  bool IsAttachedToCart() const { return attached_to_cart_; }

  void SetType(int type) { type_ = type; }
  int GetType() const { return type_; }

  Real3 GetOlderVelocity() const { return older_velocity_; }
  void SetOlderVelocity(const Real3& velocity) { older_velocity_ = velocity; }

  real_t GetOxygenConsumptionRate() const { return oxygen_consumption_rate_; }
  void SetOxygenConsumptionRate(real_t rate) {
    oxygen_consumption_rate_ = rate;
  }

  real_t GetImmunostimulatoryFactorSecretionRate() const {
    return immunostimulatory_factor_secretion_rate_;
  }
  void SetImmunostimulatoryFactorSecretionRate(real_t rate) {
    immunostimulatory_factor_secretion_rate_ = rate;
  }

  real_t GetTargetTotalVolume() const;

  /// Returns the diffusion grid for oxygen
  DiffusionGrid* GetOxygenDiffusionGrid() const { return oxygen_dgrid_; }
  /// Returns the diffusion grid for immunostimulatory factors
  DiffusionGrid* GetImmunostimulatoryFactorDiffusionGrid() const {
    return immunostimulatory_factor_dgrid_;
  }

  /// Change volume using exponential relaxation equation
  ///
  /// This method explicitly solves the system of exponential relaxation
  /// differential equations using a discrete update step. It is used to grow or
  /// shrink the volume (and proportions) smoothly toward a desired target
  /// volume over time. The relaxation rate controls the speed of convergence
  /// and dt=1 (the time_step).
  ///
  /// @param relaxation_rate_cytoplasm Relaxation rate for cytoplasm volume
  /// changes
  /// @param relaxation_rate_nucleus Relaxation rate for nucleus volume changes
  /// @param relaxation_rate_fluid Relaxation rate for fluid volume changes
  void ChangeVolumeExponentialRelaxationEquation(
      real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus,
      real_t relaxation_rate_fluid);

  /// Calculate displacement of the cell
  ///
  /// Computes the displacement of the cell based on interaction forces.
  ///
  /// @param force Pointer to the interaction force object
  /// @param squared_radius The squared radius of the cell
  /// @param dt The time step for the simulation
  /// @return The calculated displacement vector
  Real3 CalculateDisplacement(const InteractionForce* force,
                              real_t squared_radius, real_t dt) override;

  /// Consume or secrete substances
  ///
  /// Computes new oxygen or immunostimulatory factor concentration after
  /// consumption or secretion by the cell.
  ///
  /// @param substance_id The ID of the substance (oxygen or immunostimulatory
  /// factor)
  /// @param old_concentration The previous concentration of the substance
  /// @return The new concentration after consumption/secretion
  real_t ConsumeSecreteSubstance(int substance_id, real_t old_concentration);

  /// Compute constants for consumption and secretion
  ///
  /// Updates constants after the cell's change of volume or quantities.
  /// These constants are used in the consumption/secretion differential
  /// equations.
  void ComputeConstantsConsumptionSecretion();

 private:
  /// Current state of the tumor cell
  // NOLINTNEXTLINE(readability-identifier-naming)
  TumorCellState state_ = TumorCellState::kAlive;

  /// Timer to track time in the current state (in minutes)
  // NOLINTNEXTLINE(readability-identifier-naming)
  int timer_state_ = 0;

  /// Pointer to the oxygen diffusion grid
  // NOLINTNEXTLINE(readability-identifier-naming)
  DiffusionGrid* oxygen_dgrid_ = nullptr;

  /// Pointer to the immunostimulatory factor diffusion grid
  // NOLINTNEXTLINE(readability-identifier-naming)
  DiffusionGrid* immunostimulatory_factor_dgrid_ = nullptr;

  /// Level of oncoprotein expression
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t oncoproteine_level_ = 0.0;

  /// Transition random rate between states:
  /// Affects the probability of transitioning and depends on the individual
  /// cell. This rate is kept constant during the cell's lifetime.
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t transformation_random_rate_ = 0.0;

  /// Flag indicating if the cell is attached to a CAR-T cell
  // NOLINTNEXTLINE(readability-identifier-naming)
  bool attached_to_cart_ = false;

  /// Fluid fraction of the cell volume
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t fluid_fraction_ = 0.0;

  /// Volume of the nucleus
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t nuclear_volume_ = 0.0;

  /// Target cytoplasm solid volume for exponential relaxation
  ///
  /// Used for growing or shrinking tumor cells. The volume change follows
  /// an exponential relaxation equation toward this target volume.
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t target_cytoplasm_solid_ = 0.0;

  /// Target nucleus solid volume for exponential relaxation
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t target_nucleus_solid_ = 0.0;

  /// Target fluid fraction for exponential relaxation
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t target_fraction_fluid_ = 0.0;

  /// Target relation between cytoplasm and nucleus volumes
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t target_relation_cytoplasm_nucleus_ = 0.0;

  /// Cell type according to oncoprotein level:
  /// Types 1-4: 1 is the most mutated and proliferative type, 4 is the least
  /// aggressive. Type 5 means the cell is dead.
  // NOLINTNEXTLINE(readability-identifier-naming)
  int type_ = 0;

  /// Velocity of the cell in the previous time step
  // NOLINTNEXTLINE(readability-identifier-naming)
  Real3 older_velocity_ = {};

  /// Rate of oxygen consumption by the cell
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t oxygen_consumption_rate_ = 0.0;

  /// Rate of immunostimulatory factor secretion by the cell
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t immunostimulatory_factor_secretion_rate_ = 0.0;

  /// Constant 1 for oxygen consumption/secretion differential equation solution
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t constant1_oxygen_ = 0.0;

  /// Constant 2 for oxygen consumption/secretion differential equation solution
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t constant2_oxygen_ = 0.0;

  /// Constant 1 for immunostimulatory factor consumption/secretion differential
  /// equation solution
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t constant1_immunostimulatory_factor_ = 0.0;

  /// Constant 2 for immunostimulatory factor consumption/secretion differential
  /// equation solution
  // NOLINTNEXTLINE(readability-identifier-naming)
  real_t constant2_immunostimulatory_factor_ = 0.0;
};

/// Behavior class for controlling tumor cell state transitions and growth
///
/// This behavior handles the state control logic for tumor cells, managing
/// transitions between different cell states, growth, proliferation, and death
/// processes. It includes logic for determining when cells should enter
/// necrosis based on oxygen levels and other environmental factors.

struct StateControlGrowProliferate : public Behavior {
  BDM_BEHAVIOR_HEADER(StateControlGrowProliferate, Behavior, 1);

  StateControlGrowProliferate() { AlwaysCopyToNew(); }

  // Special member functions
  StateControlGrowProliferate(const StateControlGrowProliferate&) = default;
  StateControlGrowProliferate& operator=(const StateControlGrowProliferate&) = default;
  StateControlGrowProliferate(StateControlGrowProliferate&&) = default;
  StateControlGrowProliferate& operator=(StateControlGrowProliferate&&) = default;

  ~StateControlGrowProliferate() override = default;

  /// Execute the state control and growth behavior
  void Run(Agent* agent) override;

 private:
  /// Compute the probability of the cell entering necrosis
  ///
  /// Determines whether a cell should enter necrosis based on oxygen levels
  ///
  /// @param oxygen_level Current oxygen concentration at the cell's location
  /// @param cell Pointer to the tumor cell being evaluated
  /// @return True if the cell should enter necrosis, false otherwise
  bool ShouldEnterNecrosis(real_t oxygen_level, TumorCell* cell) const;
};

}  // namespace bdm

#endif  // TUMOR_CELL_H_