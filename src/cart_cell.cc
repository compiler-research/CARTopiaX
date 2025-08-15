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

#include "cart_cell.h"
#include "hyperparams.h"
#include "tumor_cell.h"
#include "utils_aux.h"
#include "core/agent/agent.h"
#include "core/agent/new_agent_event.h"
#include "core/container/math_array.h"
#include "core/diffusion/diffusion_grid.h"
#include "core/functor.h"
#include "core/interaction_force.h"
#include "core/real_t.h"
#include "core/resource_manager.h"
#include "core/util/log.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace bdm {

CarTCell::CarTCell(const Real3& position) {
  SetPosition(position);
  SetVolume(kDefaultVolumeNewCarTCell);
  const ResourceManager& rm = *Simulation::GetActive()->GetResourceManager();
  oxygen_dgrid_ = rm.GetDiffusionGrid("oxygen");
  immunostimulatory_factor_dgrid_ =
      rm.GetDiffusionGrid("immunostimulatory_factor");
  ComputeConstantsConsumptionSecretion();
}

// Cart cells can move if they are alive and not attached to a tumor cell
bool CarTCell::DoesCellMove() {
  return (state_ == CarTCellState::kAlive && !attached_to_tumor_cell_);
}

real_t CarTCell::GetTargetTotalVolume() const {
  return GetTargetNucleusSolid() * (1 + GetTargetRelationCytoplasmNucleus()) /
         (1 - GetTargetFractionFluid());
}

// This method explicitly solves the system of exponential relaxation
// differential equation using a discrete update step. It is used to shrink the
// volume (and proportions) smoothly toward a desired target volume over time
// whe the cell is apoptotic. The relaxations rate controls the speed of
// convergence
void CarTCell::ChangeVolumeExponentialRelaxationEquation(
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
Real3 CarTCell::CalculateDisplacement(const InteractionForce* force,
                                      real_t squared_radius, real_t /*dt*/) {
  // real_t h = dt;
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a
  // cleaner code instead of hardcoding it here
  squared_radius = gKSquaredMaxDistanceNeighborsForce;

  // the physics force to move the point mass

  Real3 translation_velocity_on_point_mass{0, 0, 0};

  // We check for every neighbor object if they touch us, i.e. push us
  // away and agreagate the velocities

  uint64_t non_zero_neighbor_forces = 0;
  if (!IsStatic()) {
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto calculate_neighbor_forces =
        L2F([&](Agent* neighbor, real_t /*squared_distance*/) {
          auto neighbor_force = force->Calculate(this, neighbor);
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
real_t CarTCell::ConsumeSecreteSubstance(int substance_id,
                                         real_t old_concentration) {
  real_t res = NAN;
  if (substance_id == oxygen_dgrid_->GetContinuumId()) {
    // consuming oxygen
    res = (old_concentration + constant1_oxygen_) / constant2_oxygen_;
  } else if (substance_id ==
             immunostimulatory_factor_dgrid_->GetContinuumId()) {
    // This point should never be reached
    res = old_concentration;
  } else {
    throw std::invalid_argument("Unknown substance id: " +
                                std::to_string(substance_id));
  }
  return res;
}

// Recompute Consumption constants whenever oxygen_consumption_rate_ or the
// volume changes
void CarTCell::ComputeConstantsConsumptionSecretion() {
  // constant1_= dt · (V_k / V_voxel) · S_k · ρ*_k)
  // constant2_ = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // where:
  // S_k    = secretion rate of cell k
  // U_k    = uptake (consumption) rate of cell k
  // ρ*_k   = saturation (target) density for secretion
  // V_k    = volume of the cell k
  // V_voxel = volume of the voxel containing the cell
  // dt     = simulation time step
  const real_t volume = GetVolume();
  // compute the constants for the differential equation explicit solution: for
  // oxygen and immunostimulatory factor
  // dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt
  // · (V_k / V_voxel) · S_k · ρ*_k)
  constant1_oxygen_ = 0.;
  // 1 + dt*(cell_volume/voxel_volume)*(quantity_secretion +
  // quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  //  Scale by the volume of the cell in the Voxel and time step
  constant2_oxygen_ =
      1 + kDtSubstances * (volume / kVoxelVolume) * (oxygen_consumption_rate_);
}

/// Main behavior executed at each simulation step
void StateControlCart::Run(Agent* agent) {
  auto* sim = Simulation::GetActive();
  // Run only every kDtCycle minutes, fmod does not work with the type
  // returned by GetSimulatedTime()
  if (sim->GetScheduler()->GetSimulatedSteps() % kStepsPerCycle != 0) {
    return;
  }

  if (auto* cell = dynamic_cast<CarTCell*>(agent)) {
    switch (cell->GetState()) {
      case CarTCellState::kAlive: {
        // the cell is growing to real_t its size before mitosis
        // Probability of death= 1/CurrentLiveTime, division by 0
        if (sim->GetRandom()->Uniform(1.0) <
            kDtCycle / std::max(cell->GetCurrentLiveTime(), kEpsilon)) {
          // the cell Dies
          cell->SetState(CarTCellState::kApoptotic);
          // Reset timer_state, it should be 0 anyway
          cell->SetTimerState(0);
          // Set target volume to 0 (the cell will shrink)
          cell->SetTargetCytoplasmSolid(0.0);
          cell->SetTargetNucleusSolid(0.0);
          cell->SetTargetFractionFluid(0.0);
          cell->SetTargetRelationCytoplasmNucleus(0.0);
          // Reduce oxygen consumption
          cell->SetOxygenConsumptionRate(cell->GetOxygenConsumptionRate() *
                                         kReductionConsumptionDeadCells);
          // Update constants for all Consumption of oxygen
          cell->ComputeConstantsConsumptionSecretion();
          // Detach from tumor cell if it was attached
          if (cell->IsAttachedToTumorCell()) {
            cell->GetAttachedCell()->SetAttachedToCart(false);
            cell->SetAttachedCell(nullptr);
            cell->SetAttachedToTumorCell(false);
          }
        } else {
          // decrease current life time
          cell->SetCurrentLiveTime(
              (cell->GetCurrentLiveTime() - (kDtCycle * kDtCycle)));
        }
        break;
      }
      case CarTCellState::kApoptotic: {
        cell->SetTimerState(static_cast<int>(
            static_cast<real_t>(cell->GetTimerState()) + kDtCycle));

        cell->ChangeVolumeExponentialRelaxationEquation(
            kVolumeRelaxarionRateCytoplasmApoptotic,
            kVolumeRelaxarionRateNucleusApoptotic,
            kVolumeRelaxarionRateFluidApoptotic);

        // If the timer_state exceeds the time to transition (this is a fixed
        // duration transition)
        if (kTimeApoptosis < cell->GetTimerState()) {
          // remove the cell from the simulation
          auto* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      default: {
        Log::Error("StateControlCart::Run", "Unknown CarTCellState");
        break;
      }
    }
  } else {
    Log::Error("StateControlCart::Run", "SimObject is not a CarTCell");
  }
}

}  // namespace bdm
