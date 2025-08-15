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

#ifndef CART_CELL_H_
#define CART_CELL_H_

#include "biodynamo.h"
#include "core/util/log.h"
#include "core/util/root.h"
#include "utils_aux.h"
#include "hyperparams.h"
#include "tumor_cell.h"

namespace bdm {


// ─────────────────────────────
// CartCellState Enum Definition
// ─────────────────────────────
enum class CartCellState : int {
  //living cell
  kAlive=0,

  // Apoptotic schedule of a apoptotic cell: controled death
  kApoptotic=1 // Apoptotic phase (the cell is undergoing programmed cell death characterized by cell shrinkage)
};

// ─────────────────────────────
// CartCell Class Definition
// ─────────────────────────────
class CartCell : public Cell {
  BDM_AGENT_HEADER(CartCell, Cell, 1);

  public:
  CartCell() {}
  explicit CartCell(const Real3& position);
  virtual ~CartCell() {}

  //Getters and Setters
  void SetState(CartCellState state) { state_ = state; }
  CartCellState GetState() const { return state_; }

  void SetTimerState(int timer_state) { timer_state_ = timer_state; }
  int GetTimerState() const { return timer_state_; }

  void SetFluidFraction(real_t fluid_fraction) { fluid_fraction_ = fluid_fraction; }
  real_t GetFluidFraction() const { return fluid_fraction_; }

  void SetNuclearVolume(real_t nuclear_volume) { nuclear_volume_ = nuclear_volume; }
  real_t GetNuclearVolume() const { return nuclear_volume_; }

  void SetTargetCytoplasmSolid(real_t target_cytoplasm_solid) { target_cytoplasm_solid_ = target_cytoplasm_solid; }
  real_t GetTargetCytoplasmSolid() const { return target_cytoplasm_solid_; }

  void SetTargetNucleusSolid(real_t target_nucleus_solid) { target_nucleus_solid_ = target_nucleus_solid; }
  real_t GetTargetNucleusSolid() const { return target_nucleus_solid_; }  

  void SetTargetFractionFluid(real_t target_fraction_fluid) { target_fraction_fluid_ = target_fraction_fluid; }
  real_t GetTargetFractionFluid() const { return target_fraction_fluid_; }  

  void SetTargetRelationCytoplasmNucleus(real_t target_relation_cytoplasm_nucleus) { target_relation_cytoplasm_nucleus_ = target_relation_cytoplasm_nucleus; }
  real_t GetTargetRelationCytoplasmNucleus() const { return target_relation_cytoplasm_nucleus_; }

  void SetAttachedToTumorCell(bool attached) { attached_to_tumor_cell_ = attached; }
  bool IsAttachedToTumorCell() const { return attached_to_tumor_cell_; }

  Real3 GetOlderVelocity() const { return older_velocity_; }
  void SetOlderVelocity(const Real3& velocity) { older_velocity_ = velocity; }

  real_t GetOxygenConsumptionRate() const { return oxygen_consumption_rate_; }
  void SetOxygenConsumptionRate(real_t rate) { oxygen_consumption_rate_ = rate; }

  real_t GetCurrentLiveTime() const { return current_live_time_; }
  void SetCurrentLiveTime(real_t time) { current_live_time_ = time; }

  TumorCell* GetAttachedCell() const { return attached_cell_; }
  void SetAttachedCell(TumorCell* cell) { attached_cell_ = cell; }

  //returns whether the cell moves by its own
  bool DoesCellMove();

  real_t GetTargetTotalVolume();

  /// Returns the diffusion grid for oxygen
  DiffusionGrid* GetOxygenDiffusionGrid() const { return oxygen_dgrid_; }
  /// Returns the diffusion grid for immunostimulatory factors
  DiffusionGrid* GetImmunostimulatoryFactorDiffusionGrid() const { return immunostimulatory_factor_dgrid_; }

  // This method explicitly solves the system of exponential relaxation differential equation using a discrete 
  // update step. It is used to grow or shrink the volume (and proportions) smoothly toward a desired target 
  // volume over time. The relaxations rate controls the speed of convergence and dt=1 (the time_step). 
  void ChangeVolumeExponentialRelaxationEquation(real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus, real_t relaxation_rate_fluid);

  //compute Displacement
  Real3 CalculateDisplacement(const InteractionForce* force,
                              real_t squared_radius, real_t dt) override;

  //Compute new oxygen or immunostimulatory factor concentration after consumption/ secretion
  real_t ConsumeSecreteSubstance(int substance_id, real_t old_concentration);

  //constants after cell's change of volume or quantities
  void ComputeConstantsConsumptionSecretion();

 //Attributes
 private:
  CartCellState state_;  // Current state of the cart cell
  int timer_state_;//timer to track time in the current state (in minutes): used for apoptotic state
  DiffusionGrid* oxygen_dgrid_;  // Pointer to the oxygen diffusion grid
  DiffusionGrid* immunostimulatory_factor_dgrid_;  // Pointer to the immunostimulatory_factor diffusion grid
  bool attached_to_tumor_cell_; // Flag to indicate if the cell is attached to a tumor cell
  real_t current_live_time_; // Current time untill apoptosis
  //volumes
  real_t fluid_fraction_;
  real_t nuclear_volume_; // Volume of the nucleus
  // Target volume for shrinking apoptotic cells. The change of volume follows a exponential relaxation equation with this target volume
  real_t target_cytoplasm_solid_; 
  real_t target_nucleus_solid_;
  real_t target_fraction_fluid_;
  real_t target_relation_cytoplasm_nucleus_;
  Real3 older_velocity_; // Velocity of the cell in the previous step
  real_t oxygen_consumption_rate_;
  real_t immunostimulatory_factor_secretion_rate_;
  //constants for ConsumptionSecretion differential equation solution
  real_t constant1_oxygen_;
  real_t constant2_oxygen_;
  TumorCell* attached_cell_; // Pointer to the attached tumor cell
};

// ─────────────────────────────
// Behavior: StateControlCart
// ─────────────────────────────
struct StateControlCart : public Behavior {
  BDM_BEHAVIOR_HEADER(StateControlCart, Behavior, 1);

  StateControlCart() { AlwaysCopyToNew(); }
  virtual ~StateControlCart() {}

  void Run(Agent* agent) override;
};

}  // namespace bdm

#endif  // CART_CELL_H_
