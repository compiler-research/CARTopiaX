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

#ifndef CART_CELL_H_
#define CART_CELL_H_

#include "tumor_cell.h"
#include "core/agent/agent.h"
#include "core/agent/cell.h"
#include "core/agent/new_agent_event.h"
#include "core/behavior/behavior.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/interaction_force.h"
#include "core/real_t.h"

namespace bdm {

/// Enumeration defining the possible states of a CAR-T cell
///
/// This enum class represents the different states that a CAR-T cell can be in
/// during its lifecycle in the simulation.
enum class CarTCellState : int {
  kAlive =
      0,  ///< Living cell state - the cell is alive and functioning normally

  kApoptotic =
      1  ///< Apoptotic phase - the cell is undergoing programmed cell death
         ///< characterized by cell shrinkage and controlled death
};

/// CAR-T cell class implementation
///
/// This class represents a CAR-T (Chimeric Antigen Receptor T-cell) in the
/// simulation. It inherits from the base Cell class and includes specific
/// behaviors and properties related to CAR-T cell biology, including states,
/// volume dynamics, and interactions with tumor cells.
class CarTCell : public Cell {
  // NOLINTNEXTLINE(modernize-type-traits)
  BDM_AGENT_HEADER(CarTCell, Cell, 1);

 public:
  CarTCell() = default;
  explicit CarTCell(const Real3& position);

  // Copy and move constructors/destructors/assignment operators
  CarTCell(const CarTCell&) = default;
  CarTCell(CarTCell&&) = default;
  CarTCell& operator=(const CarTCell&) = delete;
  CarTCell& operator=(CarTCell&&) = delete;
  ~CarTCell() override = default;

  /// Getters and Setters
  void SetState(CarTCellState state) { state_ = state; }
  CarTCellState GetState() const { return state_; }

  void SetTimerState(int timer_state) { timer_state_ = timer_state; }
  int GetTimerState() const { return timer_state_; }

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

  void SetAttachedToTumorCell(bool attached) {
    attached_to_tumor_cell_ = attached;
  }
  bool IsAttachedToTumorCell() const { return attached_to_tumor_cell_; }

  Real3 GetOlderVelocity() const { return older_velocity_; }
  void SetOlderVelocity(const Real3& velocity) { older_velocity_ = velocity; }

  real_t GetOxygenConsumptionRate() const { return oxygen_consumption_rate_; }
  void SetOxygenConsumptionRate(real_t rate) {
    oxygen_consumption_rate_ = rate;
  }

  real_t GetCurrentLiveTime() const { return current_live_time_; }
  void SetCurrentLiveTime(real_t time) { current_live_time_ = time; }

  TumorCell* GetAttachedCell() const { return attached_cell_; }
  void SetAttachedCell(TumorCell* cell) { attached_cell_ = cell; }

  /// Returns whether the cell moves by its own
  bool DoesCellMove();

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
  /// Current state of the CAR-T cell
  CarTCellState state_ = CarTCellState::kAlive;

  /// Timer to track time in the current state (in minutes)
  /// Used for apoptotic state timing
  int timer_state_ = 0;

  /// Pointer to the oxygen diffusion grid
  DiffusionGrid* oxygen_dgrid_ = nullptr;

  /// Pointer to the immunostimulatory factor diffusion grid
  DiffusionGrid* immunostimulatory_factor_dgrid_ = nullptr;

  /// Flag indicating if the cell is attached to a tumor cell
  bool attached_to_tumor_cell_ = false;

  /// Current time until apoptosis
  real_t current_live_time_ = 0.0;

  /// Fluid fraction of the cell volume
  real_t fluid_fraction_ = 0.0;

  /// Volume of the nucleus
  real_t nuclear_volume_ = 0.0;

  /// Target cytoplasm solid volume for exponential relaxation
  /// Used during volume changes following exponential relaxation equation
  real_t target_cytoplasm_solid_ = 0.0;

  /// Target nucleus solid volume for exponential relaxation
  real_t target_nucleus_solid_ = 0.0;

  /// Target fluid fraction for exponential relaxation
  real_t target_fraction_fluid_ = 0.0;

  /// Target relation between cytoplasm and nucleus volumes
  real_t target_relation_cytoplasm_nucleus_ = 0.0;

  /// Velocity of the cell in the previous time step
  Real3 older_velocity_ = {0.0, 0.0, 0.0};

  /// Rate of oxygen consumption by the cell
  real_t oxygen_consumption_rate_ = 0.0;

  /// Rate of immunostimulatory factor secretion by the cell
  real_t immunostimulatory_factor_secretion_rate_ = 0.0;

  /// Constant 1 for oxygen consumption/secretion differential equation solution
  real_t constant1_oxygen_ = 0.0;

  /// Constant 2 for oxygen consumption/secretion differential equation solution
  real_t constant2_oxygen_ = 0.0;

  /// Pointer to the attached tumor cell
  TumorCell* attached_cell_ = nullptr;
};

/// Behavior class for controlling CAR-T cell state transitions
///
/// This behavior handles the state control logic for CAR-T cells, managing
/// transitions between different cell states: alive and apoptotic phases.
/// It inherits from the base Behavior class and implements the Run method to
/// execute the state control logic during simulation steps.
struct StateControlCart : public Behavior {
  // NOLINTNEXTLINE(llvm-else-after-return, moderinize-type-traits,
  // cppcoreguidelines-owning-memory)
  BDM_BEHAVIOR_HEADER(StateControlCart, Behavior, 1);

  StateControlCart() { AlwaysCopyToNew(); }

  // Special member functions
  StateControlCart(const StateControlCart&) = default;
  StateControlCart(StateControlCart&&) = default;
  StateControlCart& operator=(const StateControlCart&) = default;
  StateControlCart& operator=(StateControlCart&&) = default;
  ~StateControlCart() override = default;

  /// Execute the state control behavior
  void Run(Agent* agent) override;
};

}  // namespace bdm

#endif  // CART_CELL_H_
