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

#include "tumor_cell.h"
#include "hyperparams.h"
#include "utils_aux.h"
#include "core/agent/agent.h"
#include "core/agent/cell_division_event.h"
#include "core/agent/new_agent_event.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/functor.h"
#include "core/interaction_force.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include "core/simulation.h"
#include "core/util/log.h"
#include "core/util/random.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace bdm {

TumorCell::TumorCell(const Real3& position) {
  SetPosition(position);

  // volumes
  // Set default volume
  SetVolume(kDefaultVolumeNewTumorCell);
  // Set default fluid fraction
  SetFluidFraction(kDefaultFractionFluidTumorCell);
  // Set default nuclear volume
  SetNuclearVolume(kDefaultVolumeNucleusTumorCell);
  // target volumes
  // Set target fraction of fluid
  SetTargetFractionFluid(kDefaultFractionFluidTumorCell);
  // Set target relation between cytoplasm and nucleus
  SetTargetRelationCytoplasmNucleus(
      (kDefaultVolumeNewTumorCell - kDefaultVolumeNucleusTumorCell) /
      (kEpsilon + kDefaultVolumeNucleusTumorCell));
  // Set target nucleus solid volume
  SetTargetNucleusSolid(kDefaultVolumeNucleusTumorCell *
                        (1 - kDefaultFractionFluidTumorCell));
  // Set target cytoplasm solid volume
  SetTargetCytoplasmSolid(
      (kDefaultVolumeNewTumorCell - kDefaultVolumeNucleusTumorCell) *
      (1 - kDefaultFractionFluidTumorCell));

  // Set initial oncoprotein level with a truncated normal distribution
  SetOncoproteinLevel(
      SamplePositiveGaussian(kOncoproteinMean, kOncoproteinStandardDeviation));
  ResourceManager* rm = Simulation::GetActive()->GetResourceManager();
  // Pointer to oxygen diffusion grid
  oxygen_dgrid_ = rm->GetDiffusionGrid("oxygen");
  // Pointer to immunostimulatory_factor diffusion grid
  immunostimulatory_factor_dgrid_ =
      rm->GetDiffusionGrid("immunostimulatory_factor");
  // Set state transition random rate
  SetTransformationRandomRate();

  // Add Consumption and Secretion
  // Set default oxygen consumption rate
  SetOxygenConsumptionRate(kDefaultOxygenConsumption);
  // Set default immunostimulatory factor secretion rate
  SetImmunostimulatoryFactorSecretionRate(
      kRateSecretionImmunostimulatoryFactor);
  // Compute constants for all ConsumptionSecretion of Oxygen and
  // Immunostimulatory Factor
  ComputeConstantsConsumptionSecretion();
}

/// Called when a new agent is created (e.g., after cell division)
void TumorCell::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);

  // if the cell is created from division
  if (auto* mother = dynamic_cast<TumorCell*>(event.existing_agent)) {
    if (event.GetUid() == CellDivisionEvent::kUid) {
      // Initialize daughter cell from mother cell
      // state after division
      state_ = TumorCellState::kAlive;
      timer_state_ = 0;
      // diffusion grids
      // Pointer to the oxygen diffusion grid
      oxygen_dgrid_ = mother->oxygen_dgrid_;
      // Pointer to the immunostimulatory_factor diffusion grid
      immunostimulatory_factor_dgrid_ = mother->immunostimulatory_factor_dgrid_;
      // inherit oncoprotein level from mother cell
      this->SetOncoproteinLevel(mother->oncoprotein_level_);
      // inherit oxygen consumption from mother cell
      this->SetOxygenConsumptionRate(mother->GetOxygenConsumptionRate());
      // inherit immunostimulatory factor secretion rate from mother cell
      this->SetImmunostimulatoryFactorSecretionRate(
          mother->GetImmunostimulatoryFactorSecretionRate());

      // Update the constants for all ConsumptionSecretion
      mother->ComputeConstantsConsumptionSecretion();
      this->ComputeConstantsConsumptionSecretion();

      // divide mother's nuclear volume by 2
      const real_t new_nuclear_volume = mother->GetNuclearVolume() / kHalf;
      // Set mother's nuclear volume to the new volume
      mother->SetNuclearVolume(new_nuclear_volume);
      this->SetNuclearVolume(new_nuclear_volume);

      // Inherit mother's fluid fraction and velocity
      // Set fluid fraction to mother's fluid fraction
      this->SetFluidFraction(mother->GetFluidFraction());
      // Copy velocity from mother cell
      this->SetOlderVelocity(mother->GetOlderVelocity());

      // inherit target volumes of the daughter cell
      this->SetTargetFractionFluid(mother->GetTargetFractionFluid());
      this->SetTargetRelationCytoplasmNucleus(
          mother->GetTargetRelationCytoplasmNucleus());
      this->SetTargetNucleusSolid(mother->GetTargetNucleusSolid());
      this->SetTargetCytoplasmSolid(mother->GetTargetCytoplasmSolid());

      // Set state transition random rate
      this->SetTransformationRandomRate();
      // Initially not attached to a cart
      this->attached_to_cart_ = false;
    }
  }
}

void TumorCell::SetOncoproteinLevel(real_t level) {
  oncoprotein_level_ = level;
  // cell type
  if (level >= kThresholdCancerCellType1) {
    // between 1.5 and 2.0
    type_ = TumorCellType::kType1;
  } else if (level >= kThresholdCancerCellType2 &&
             level < kThresholdCancerCellType1) {
    type_ = TumorCellType::kType2;
  } else if (level >= kThresholdCancerCellType3 &&
             level < kThresholdCancerCellType2) {
    type_ = TumorCellType::kType3;
  } else if (level >= kThresholdCancerCellType4 &&
             level < kThresholdCancerCellType3) {
    type_ = TumorCellType::kType4;
  } else {
    // undefined type
    type_ = TumorCellType::kType0;
  }
}

void TumorCell::SetTransformationRandomRate() {
  // avoid division by zero with kEpsilon
  transformation_random_rate_ =
      1 / (std::max(SamplePositiveGaussian(
                        kAverageTimeTransformationRandomRate,
                        kStandardDeviationTransformationRandomRate) *
                        kMinutesInAnHour,
                    kEpsilon));
}

real_t TumorCell::GetTargetTotalVolume() const {
  return GetTargetNucleusSolid() * (1 + GetTargetRelationCytoplasmNucleus()) /
         (1 - GetTargetFractionFluid());
}

// This method explicitly solves the system of exponential relaxation
// differential equation using a discrete update step. It is used to grow or
// shrink the volume (and proportions) smoothly toward a desired target volume
// over time. The relaxations rate controls the speed of convergence
void TumorCell::ChangeVolumeExponentialRelaxationEquation(
    real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus,
    real_t relaxation_rate_fluid) {
  // Exponential relaxation towards the target volume
  const real_t current_total_volume = GetVolume();
  const real_t fluid_fraction = GetFluidFraction();
  const real_t nuclear_volume = GetNuclearVolume();

  const real_t current_nuclear_solid = nuclear_volume * (1 - fluid_fraction);
  const real_t current_cytoplasm_solid =
      (current_total_volume - nuclear_volume) * (1 - fluid_fraction);

  const real_t current_fluid = fluid_fraction * current_total_volume;

  // Update fluid volume
  real_t new_fluid =
      current_fluid +
      kDtCycle * relaxation_rate_fluid *
          (GetTargetFractionFluid() * current_total_volume - current_fluid);
  // Clamp to zero to prevent negative volumes
  if (new_fluid < 0.0) {
    new_fluid = 0.0;
  }

  const real_t nuclear_fluid =
      new_fluid * (nuclear_volume / current_total_volume);
  // real_t cytoplasm_fluid = new_fluid - nuclear_fluid;

  real_t nuclear_solid = current_nuclear_solid +
                         kDtCycle * relaxation_rate_nucleus *
                             (GetTargetNucleusSolid() - current_nuclear_solid);
  // Clamp to zero to prevent negative volumes
  if (nuclear_solid < 0.0) {
    nuclear_solid = 0.0;
  }

  const real_t target_cytoplasm_solid =
      GetTargetRelationCytoplasmNucleus() * GetTargetNucleusSolid();
  real_t cytoplasm_solid =
      current_cytoplasm_solid +
      kDtCycle * relaxation_rate_cytoplasm *
          (target_cytoplasm_solid - current_cytoplasm_solid);
  // Clamp to zero to prevent negative volumes
  if (cytoplasm_solid < 0.0) {
    cytoplasm_solid = 0.0;
  }

  const real_t new_total_solid = nuclear_solid + cytoplasm_solid;

  const real_t total_nuclear = nuclear_solid + nuclear_fluid;

  // real_t total_cytoplasm= cytoplasm_solid + cytoplasm_fluid;

  const real_t new_volume = new_total_solid + new_fluid;

  // Avoid division by zero
  const real_t new_fraction_fluid = new_fluid / (kEpsilon + new_volume);

  // Update the cell's properties
  // if the volume has changed
  if (new_volume != current_total_volume) {
    SetVolume(new_volume);
    // Update constants for all ConsumptionSecretion of Oxygen and
    // Immunostimulatory Factors
    ComputeConstantsConsumptionSecretion();
  }
  SetFluidFraction(new_fraction_fluid);
  SetNuclearVolume(total_nuclear);
}

// compute Displacement
Real3 TumorCell::CalculateDisplacement(const InteractionForce* force,
                                       real_t squared_radius, real_t /*dt*/) {
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a
  // cleaner code instead of hardcoding it here
  squared_radius = gKSquaredMaxDistanceNeighborsForce;

  Real3 translation_velocity_on_point_mass{0, 0, 0};

  // We check for every neighbor object if they touch us, i.e. push us
  // away and agreagate the velocities

  uint64_t non_zero_neighbor_forces = 0;
  if (!IsStatic()) {
    ExecutionContext* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto calculate_neighbor_forces =
        L2F([&](Agent* neighbor, real_t /*squared_distance*/) {
          Real4 neighbor_force = force->Calculate(this, neighbor);
          if (neighbor_force[0] != 0 || neighbor_force[1] != 0 ||
              neighbor_force[2] != 0) {
            non_zero_neighbor_forces++;
            translation_velocity_on_point_mass[0] += neighbor_force[0];
            translation_velocity_on_point_mass[1] += neighbor_force[1];
            translation_velocity_on_point_mass[2] += neighbor_force[2];
          }
        });
    ctxt->ForEachNeighbor(calculate_neighbor_forces, *this, squared_radius);

    if (non_zero_neighbor_forces > 1) {
      SetStaticnessNextTimestep(false);
    }
  }

  // Two step Adams-Bashforth approximation of the time derivative for position
  // position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 *
  // velocity(t - dt) ]
  movement_at_next_step +=
      translation_velocity_on_point_mass * kDnew + older_velocity_ * kDold;

  older_velocity_ = translation_velocity_on_point_mass;

  // Displacement
  return movement_at_next_step;
}

// Compute new oxygen or immunostimulatory factor concentration after
// consumption/ secretion
real_t TumorCell::ConsumeSecreteSubstance(int substance_id,
                                          real_t old_concentration) {
  real_t res = 0.0;
  if (substance_id == oxygen_dgrid_->GetContinuumId()) {
    // consuming oxygen
    res = (old_concentration + constant1_oxygen_) / constant2_oxygen_;
  } else if (substance_id ==
             immunostimulatory_factor_dgrid_->GetContinuumId()) {
    // secreting immunostimulatory factor
    res = (old_concentration + constant1_immunostimulatory_factor_) /
          constant2_immunostimulatory_factor_;
  } else {
    throw std::invalid_argument("Unknown substance id: " +
                                std::to_string(substance_id));
  }
  return res;
}

void TumorCell::ComputeConstantsConsumptionSecretion() {
  // constant1_= dt · (V_k / V_voxel) · S_k · ρ*_k)
  // constant2_ = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // where:
  // S_k    = secretion rate of cell k
  // U_k    = uptake (consumption) rate of cell k
  // ρ*_k   = saturation (target) density for secretion
  // V_k    = volume of the cell k
  // V_voxel = volume of the voxel containing the cell
  // dt     = simulation time step

  const real_t new_volume = GetVolume();
  // compute the constants for the differential equation explicit solution: for
  // oxygen and immunostimulatory factor
  // dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt
  // · (V_k / V_voxel) · S_k · ρ*_k)
  constant1_oxygen_ = 0.;
  // Scale by the volume of the cell in the Voxel and time step
  constant1_immunostimulatory_factor_ =
      immunostimulatory_factor_secretion_rate_ *
      kSaturationDensityImmunostimulatoryFactor * kDtSubstances *
      (new_volume / kVoxelVolume);
  // 1 + dt*(cell_volume/voxel_volume)*(quantity_secretion +
  // quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  //  Scale by the volume of the cell in the Voxel and time step
  constant2_oxygen_ = 1 + kDtSubstances * (new_volume / kVoxelVolume) *
                              (oxygen_consumption_rate_);
  // Scale by the volume of the cell in the Voxel and time step
  constant2_immunostimulatory_factor_ =
      1 + kDtSubstances * (new_volume / kVoxelVolume) *
              (immunostimulatory_factor_secretion_rate_);
}

void TumorCell::StartApoptosis() {
  // If the cell is already dead, do nothing
  if (IsDead()) {
    return;
  }

  // The cell Dies
  SetState(TumorCellState::kApoptotic);

  // Reset timer_state
  SetTimerState(0);
  // Set type to 5 to indicate dead cell
  SetType(TumorCellType::kType5);
  // Set target volume to 0 (the cell shrinks)
  SetTargetCytoplasmSolid(0.0);
  SetTargetNucleusSolid(0.0);
  SetTargetFractionFluid(0.0);
  SetTargetRelationCytoplasmNucleus(0.0);
  // Reduce oxygen consumption
  SetOxygenConsumptionRate(GetOxygenConsumptionRate() *
                           kReductionConsumptionDeadCells);
  // Stop Immunostimulatory Factor Secretion
  SetImmunostimulatoryFactorSecretionRate(0.0);
  // Update constants for consumption/secretion differential equation solving
  ComputeConstantsConsumptionSecretion();
}

/// Main behavior executed at each simulation step
void StateControlGrowProliferate::Run(Agent* agent) {
  Simulation* sim = Simulation::GetActive();
  // Run only every kDtCycle minutes, fmod does not work with the type
  // returned by GetSimulatedTime()
  if (sim->GetScheduler()->GetSimulatedSteps() % kStepsPerCycle != 0) {
    return;
  }

  if (auto* cell = dynamic_cast<TumorCell*>(agent)) {
    if (cell->IsAttachedToCart()) {
      // If the cell is attached to a cart, skip the state control and growth
      return;
    }
    // Oxygen levels
    const Real3 current_position = cell->GetPosition();
    // Pointer to the oxygen diffusion grid
    DiffusionGrid* oxygen_dgrid = cell->GetOxygenDiffusionGrid();
    const real_t oxygen_level = oxygen_dgrid->GetValue(current_position);

    switch (cell->GetState()) {
      case TumorCellState::kAlive: {
        // the cell is growing to real_t its size before mitosis
        // Increase timer_state to track time in this state (kDtCycle minutes
        // per step)
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);

        // Enter necrosis if oxygen level is too low
        if (ShouldEnterNecrosis(oxygen_level, cell)) {
          // Exit the function to prevent further processing
          return;
        }
        ManageLivingCell(cell, oxygen_level);
        break;
      }
      case TumorCellState::kNecroticSwelling: {
        // the cell is swelling before lysing
        // Increase timer_state to track time in this state (kDtCycle minutes
        // per step)
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);
        // volume change
        //  The cell swells
        cell->ChangeVolumeExponentialRelaxationEquation(
            kVolumeRelaxarionRateCytoplasmNecroticSwelling,
            kVolumeRelaxarionRateNucleusNecroticSwelling,
            kVolumeRelaxarionRateFluidNecroticSwelling);

        // If the cell has swollen to 2 times its original volume, it lyses
        if (cell->GetVolume() >= 2 * kDefaultVolumeNewTumorCell) {
          // Change state to necrotic lysed
          cell->SetState(TumorCellState::kNecroticLysed);
          // Reset timer_state
          cell->SetTimerState(0);
          // Set target volume to 0 (the cell will shrink)
          cell->SetTargetCytoplasmSolid(0.0);
          cell->SetTargetNucleusSolid(0.0);
          cell->SetTargetFractionFluid(0.0);
          cell->SetTargetRelationCytoplasmNucleus(0.0);
          // Stop secretion and consumption rate
          // Stop consumption
          cell->SetOxygenConsumptionRate(0.0);
          // Stop secretion
          cell->SetImmunostimulatoryFactorSecretionRate(0.0);
          // Update constants for all ConsumptionSecretion of Oxygen and
          // Immunostimulatory Factors
          cell->ComputeConstantsConsumptionSecretion();
        }
        break;
      }
      case TumorCellState::kNecroticLysed: {
        // the cell is shrinking and will be removed after a certain time
        // Increase timer_state to track time in this state (kDtCycle minutes
        // per step)
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);
        // volume change
        //  The cell shrinks
        cell->ChangeVolumeExponentialRelaxationEquation(
            kVolumeRelaxarionRateCytoplasmNecroticLysed,
            kVolumeRelaxarionRateNucleusNecroticLysed,
            kVolumeRelaxarionRateFluidNecroticLysed);
        // If the timer_state exceeds the time to transition (this is a fixed
        // duration transition)
        if (kTimeLysis < cell->GetTimerState()) {
          // remove the cell from the simulation
          ExecutionContext* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      case TumorCellState::kApoptotic: {
        // Increase timer_state to track time in this state (kDtCycle minutes
        // per step)
        cell->SetTimerState(cell->GetTimerState() + kDtCycle);

        //  The cell shrinks
        cell->ChangeVolumeExponentialRelaxationEquation(
            kVolumeRelaxarionRateCytoplasmApoptotic,
            kVolumeRelaxarionRateNucleusApoptotic,
            kVolumeRelaxarionRateFluidApoptotic);
        // If the timer_state exceeds the time to transition (this is a fixed
        // duration transition)
        if (kTimeApoptosis < cell->GetTimerState()) {
          // remove the cell from the simulation
          ExecutionContext* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      default: {
        Log::Error("StateControlGrowProliferate::Run",
                   "Unknown TumorCellState");
        break;
      }
    }
  } else {
    Log::Error("StateControlGrowProliferate::Run",
               "SimObject is not a TumorCell");
  }
}

// ManageLivingCell function to handle living cell behavior
void StateControlGrowProliferate::ManageLivingCell(TumorCell* cell,
                                                   real_t oxygen_level) {
  // Initialize multiplier
  real_t multiplier = 1.0;
  // volume change
  cell->ChangeVolumeExponentialRelaxationEquation(
      kVolumeRelaxarionRateAliveCytoplasm, kVolumeRelaxarionRateAliveNucleus,
      kVolumeRelaxarionRateAliveFluid);
  // cell state control
  // oxygen threshold for considering an effect on the proliferation cycle
  if (oxygen_level < kOxygenSaturationInProliferation) {
    multiplier =
        (oxygen_level - kOxygenLimitForProliferation) /
        (kOxygenSaturationInProliferation - kOxygenLimitForProliferation);
  }
  // If oxygen is below the limit, set multiplier to 0
  if (oxygen_level < kOxygenLimitForProliferation) {
    multiplier = 0.0;
  }
  // Calculate the rate of state change based on oxygen level and oncoprotein
  // (min^-1)
  const real_t final_rate_transition = cell->GetTransformationRandomRate() *
                                       multiplier * cell->GetOncoproteinLevel();

  // Calculate the time to wait (in minutes)
  real_t time_to_wait = kTimeTooLarge;
  if (final_rate_transition > 0) {
    // Calculate the time to transition (in minutes)
    time_to_wait = 1. / final_rate_transition;
  }
  // If the timer_state exceeds the time to transition (this is a fixed duration
  // transition)
  if (time_to_wait < cell->GetTimerState()) {
    // mitosis: cell divides
    cell->SetState(TumorCellState::kAlive);
    cell->Divide();
    // Reset timer_state
    cell->SetTimerState(0);
  }
}

// computes the probability of the cell entering necrosis
bool StateControlGrowProliferate::ShouldEnterNecrosis(real_t oxygen_level,
                                                      TumorCell* cell) {
  // necrosis probability
  // Default multiplier for necrosis probability
  real_t multiplier = 0.0;
  // oxygen threshold for considering necrosis
  if (oxygen_level < kOxygenLimitForNecrosis) {
    multiplier = (kOxygenLimitForNecrosis - oxygen_level) /
                 (kOxygenLimitForNecrosis - kOxygenLimitForNecrosisMaximum);
  }
  // threshold for maximum necrosis probability
  if (oxygen_level < kOxygenLimitForNecrosisMaximum) {
    multiplier = 1.0;
  }
  // Calculate the probability of necrosis based on oxygen level
  // and multiply by kDtCycle since each timestep is kDtCycle minutes
  const real_t probability_necrosis =
      kDtCycle * kMaximumNecrosisRate * multiplier;

  Simulation* sim = Simulation::GetActive();
  Random* random = sim->GetRandom();
  const bool enter_necrosis = random->Uniform(0, 1) < probability_necrosis;
  // If the random number is less than the probability, enter necrosis
  if (enter_necrosis) {
    // If oxygen is too low, enter necrosis
    cell->SetState(TumorCellState::kNecroticSwelling);
    // Reset timer_state
    cell->SetTimerState(0);

    // Stop Secretion and reduce consumption
    //  Stop secretion
    cell->SetImmunostimulatoryFactorSecretionRate(0.0);
    // Reduce consumption
    cell->SetOxygenConsumptionRate(cell->GetOxygenConsumptionRate() *
                                   kReductionConsumptionDeadCells);
    // Update constants for all ConsumptionSecretion of Oxygen and
    // Immunostimulatory Factors
    cell->ComputeConstantsConsumptionSecretion();

    // The cell will swell getting filled with fluid
    cell->SetTargetCytoplasmSolid(0);
    cell->SetTargetNucleusSolid(0);
    // Set target fraction of fluid to 1.0
    cell->SetTargetFractionFluid(1.0);
    cell->SetTargetRelationCytoplasmNucleus(0.0);
    // Set type to 5 to indicate dead cell
    cell->SetType(TumorCellType::kType5);
  }
  return enter_necrosis;  // Return whether the cell entered necrosis
}

}  // namespace bdm
