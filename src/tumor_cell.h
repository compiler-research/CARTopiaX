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

#ifndef TUMOR_CELL_H_
#define TUMOR_CELL_H_

#include "biodynamo.h"
#include "core/util/log.h"
#include "core/util/root.h"
#include "hyperparams.h"
#include "utils_aux.h"

namespace bdm {

// ─────────────────────────────
// TumorCellState Enum Definition
// ─────────────────────────────
/// Enum representing the different states of a tumor cell including the different proliferation phases of the protein Ki67 expression
enum class TumorCellState : int {
  // growing schedule of a living cell
  kAlive=0,

  // Death schedule of a necrotic cell: death by necrosis
  kNecroticSwelling = 1,  // Necrotic swelling phase (the cell loses membrane integrity and starts absorbing fluid, swelling abnormally in volume before rupture)
  kNecroticLysed = 2,  // Necrotic lysed phase (the cell membrane breaks apart, releasing its contents; the cell is now considered dead and will be removed from the simulation after a defined time)
    
  // Apoptotic schedule of a apoptotic cell: controled death
  kApoptotic=3 // Apoptotic phase (the cell is undergoing programmed cell death characterized by cell shrinkage)
};

// ─────────────────────────────
// TumorCell Class Definition
// ─────────────────────────────
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);

  public:
  TumorCell() {}
  explicit TumorCell(const Real3& position);
  virtual ~TumorCell() {}

  /// Called when a new agent is created (e.g., after cell division)
  void Initialize(const NewAgentEvent& event) override;

  //Getters and Setters
  void SetState(TumorCellState state) { state_ = state; }
  TumorCellState GetState() const { return state_; }

  void SetTimerState(int timer_state) { timer_state_ = timer_state; }
  int GetTimerState() const { return timer_state_; }

  void SetOncoproteineLevel(real_t level);
  real_t GetOncoproteineLevel() const { return oncoproteine_level_; }  

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

  void SetTransformationRandomRate();
  real_t GetTransformationRandomRate() const { return transformation_random_rate_; }

  void SetAttachedToCart(bool attached) { attached_to_cart_ = attached; }
  bool IsAttachedToCart() const { return attached_to_cart_; }

  void SetType(int type) { type_ = type; }
  int GetType() const { return type_; }

  Real3 GetOlderVelocity() const { return older_velocity_; }
  void SetOlderVelocity(const Real3& velocity) { older_velocity_ = velocity; }

  real_t GetOxygenConsumptionRate() const { return oxygen_consumption_rate_; }
  void SetOxygenConsumptionRate(real_t rate) { oxygen_consumption_rate_ = rate; }

  real_t GetImmunostimulatoryFactorSecretionRate() const { return immunostimulatory_factor_secretion_rate_; }
  void SetImmunostimulatoryFactorSecretionRate(real_t rate) { immunostimulatory_factor_secretion_rate_ = rate; }

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
  TumorCellState state_;  // Current state of the tumor cell
  int timer_state_;//timer to track time in the current state (in minutes)
  DiffusionGrid* oxygen_dgrid_;  // Pointer to the oxygen diffusion grid
  DiffusionGrid* immunostimulatory_factor_dgrid_;  // Pointer to the immunostimulatory_factor diffusion grid
  real_t oncoproteine_level_;// Level of oncoproteine expression
  real_t transformation_random_rate_; // Transition random rate between states. Affects the probability of transitioning and depends on the cell: it is kept constant during the cell's life 
  bool attached_to_cart_; // Flag to indicate if the cell is attached to a cart
  //volumes
  real_t fluid_fraction_;
  real_t nuclear_volume_; // Volume of the nucleus
  // Target volume for growing (or shrinking) tumor cells. The change of volume follows a exponential relaxation equation with this target volume
  real_t target_cytoplasm_solid_; 
  real_t target_nucleus_solid_;
  real_t target_fraction_fluid_;
  real_t target_relation_cytoplasm_nucleus_;
  int type_;//type acording to the oncoproteine level: 1, 2, 3 or 4. 1 is the most muttated and ploriferative type and 4 is the least aggressive one. Type 5 means dead
  Real3 older_velocity_; // Velocity of the cell in the previous step
  real_t oxygen_consumption_rate_;
  real_t immunostimulatory_factor_secretion_rate_;
  //constants for ConsumptionSecretion differential equation solution
  real_t constant1_oxygen_;
  real_t constant2_oxygen_;
  real_t constant1_immunostimulatory_factor_;
  real_t constant2_immunostimulatory_factor_;
};

// ─────────────────────────────
// Behavior: StateControlGrowProliferate
// ─────────────────────────────
struct StateControlGrowProliferate : public Behavior {
  BDM_BEHAVIOR_HEADER(StateControlGrowProliferate, Behavior, 1);

  StateControlGrowProliferate() { AlwaysCopyToNew(); }
  virtual ~StateControlGrowProliferate() {}

  void Run(Agent* agent) override;

  private:
  // computes the probability of the cell entering necrosis
  bool ShouldEnterNecrosis(real_t oxygen_level,TumorCell* cell) const;
};

}  // namespace bdm

#endif  // TUMOR_CELL_H_