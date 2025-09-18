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
#include "core/agent/agent_pointer.h"
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

  SetCurrentLiveTime(kAverageMaximumTimeUntillApoptosisCart);
  // Add Consumption and Secretion
  //  Set default oxygen consumption rate
  SetOxygenConsumptionRate(kDefaultOxygenConsumptionCarT);
  // Compute constants for all ConsumptionSecretion of Oxygen
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
  Simulation* sim = Simulation::GetActive();

  // real_t h = dt;
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a
  // cleaner code instead of hardcoding it here
  squared_radius = gKSquaredMaxDistanceNeighborsForce;

  // the physics force to move the point mass

  Real3 translation_velocity_on_point_mass{0, 0, 0};

  //--------------------------------------------
  // CAR-T self motility (in case of migration)
  //--------------------------------------------
  Real3 current_position = GetPosition();
  ExecutionContext* ctxt = sim->GetExecutionContext();
  Random* rng = sim->GetRandom();
  Real3 motility;
  if (DoesCellMove()) {
    // compute motility
    if (rng->Uniform(0.0, 1.0) < kMotilityProbability) {
      // random direction as unitary vector
      const Real3 random_direction = GenerateRandomDirection();
      Real3 direction_to_immunostimulatory_factor;
      // returns normalized gradient towards the immunostimulatory factor source
      immunostimulatory_factor_dgrid_->GetGradient(
          current_position, &direction_to_immunostimulatory_factor, true);
      // motility = bias * direction_to_immunostimulatory_factor +
      // (1-bias)*random_direction
      motility = kMigrationBiasCart * direction_to_immunostimulatory_factor +
                 kMigrationOneMinusBiasCart * random_direction;
      const real_t motility_norm_squared = motility[0] * motility[0] +
                                           motility[1] * motility[1] +
                                           motility[2] * motility[2];
      // Convert to unit direction
      if (motility_norm_squared > 0) {
        motility.Normalize();
      }
      // Scale by migration speed and add to the velocity
      translation_velocity_on_point_mass += motility * kMigrationSpeedCart;
    }
  }

  //--------------------------------------------
  // CAR-T killing or victim cell escaping
  //--------------------------------------------
  // If cell is not apoptotic
  if (state_ == CarTCellState::kAlive) {
    // if it is attached to tumor cell
    if (attached_to_tumor_cell_) {
      // try to kill the cancer cell and in case of failure see if it manages to
      // scape the order needs to be this one because it should try to kill
      // before seeing if it scapes
      if (TryToInduceApoptosis(attached_cell_ptr_, rng) ||
          rng->Uniform(0.0, 1.0) < kProbabilityEscape) {
        // the cancer cell is detached
        attached_cell_ptr_->SetAttachedToCart(false);
        // empty ID
        attached_cell_ptr_ = nullptr;
        attached_to_tumor_cell_ = false;
      }
    }

    //--------------------------------------------
    // CAR-T adhesion to victim cell
    //--------------------------------------------
    // Compute forces between the cells and check for a new attachment
    auto calculate_forces_and_elastic_displacement =
        L2F([&](Agent* neighbor, real_t /*squared_distance*/) {
          // Adhesion repulsion forces between cells
          //  We check for every neighbor object if they touch us, i.e. push us
          //  away and aggregate the velocities
          Real4 neighbor_force = force->Calculate(this, neighbor);
          translation_velocity_on_point_mass[0] += neighbor_force[0];
          translation_velocity_on_point_mass[1] += neighbor_force[1];
          translation_velocity_on_point_mass[2] += neighbor_force[2];

          // CAR-T adhesion to new victim cell
          Real3 displac = neighbor->GetPosition() - current_position;

          if (auto* cancer_cell = dynamic_cast<TumorCell*>(neighbor)) {
            // movement towards the tumor cells
            const real_t sq_norm_displac = displac[0] * displac[0] +
                                           displac[1] * displac[1] +
                                           displac[2] * displac[2];

            // The cart moves towards the tumor cell only if they are not
            // touching already If they are too close the only force affecting
            // is the adhesion force to avoid CAR-T non-stop pushing tumor
            // cells. In case of being closer than
            // gKMaxSquaredDistanceCartMovingTowardsTumorCell there is a
            // probability kProbabilityPushing for the CAR-T to keep pushing the
            // tumor cell
            if (sq_norm_displac >
                gKMaxSquaredDistanceCartMovingTowardsTumorCell) {
              translation_velocity_on_point_mass[0] +=
                  displac[0] * kElasticConstantCart;
              translation_velocity_on_point_mass[1] +=
                  displac[1] * kElasticConstantCart;
              translation_velocity_on_point_mass[2] +=
                  displac[2] * kElasticConstantCart;
            }

            // If the CAR-T has not succeeded in attaching to a tumor cell yet,
            // it tries again
            if (!attached_to_tumor_cell_) {
              TryToGetAttachedTo(cancer_cell, sq_norm_displac, rng);
            }
          }
        });
    ctxt->ForEachNeighbor(calculate_forces_and_elastic_displacement, *this,
                          squared_radius);
  }

  //--------------------------------------------
  // Two step Adams-Bashforth approximation of the time derivative for position
  // position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 *
  // velocity(t - dt) ]
  //--------------------------------------------
  movement_at_next_step +=
      translation_velocity_on_point_mass * kDnew + older_velocity_ * kDold;

  older_velocity_ = translation_velocity_on_point_mass;

  // Displacement
  return movement_at_next_step;
}

// Try to get attached to a tumor cell
void CarTCell::TryToGetAttachedTo(TumorCell* victim, real_t squared_distance,
                                  Random* rng) {
  // If the tumor cell is not already attached to a CAR-T cell, is not dead and
  // is not too far away.
  if (!victim->IsAttachedToCart() && !victim->IsDead() &&
      squared_distance < kSquaredMaxAdhesionDistanceCart) {
    // factor of how high is the oncoprotein level of the cancer cell
    real_t oncoprotein_scale_factor =
        (victim->GetOncoproteinLevel() - kOncoproteinLimit) /
        kOncoproteinDifference;
    // Clamp scale_factor to be in [0,1]
    if (oncoprotein_scale_factor > 1.0) {
      oncoprotein_scale_factor = 1.0;
    }
    // If oncoprotein level is lower than the limit the cancer cell does not get
    // detected
    if (oncoprotein_scale_factor <= 0.0) {
      // oncoprotein_scale_factor = 0.0; the probability is going to be 0 so
      // return the function is the most efficient
      return;
    }

    // factor of how far the cancer cell is
    real_t distance_scale_factor =
        (kMaxAdhesionDistanceCart - std::sqrt(squared_distance)) /
        kDifferenceCartAdhesionDistances;
    // Clamp scale_factor to be in [0,1]. We already checked that it is > 0
    // because squared_distance < kSquaredMaxAdhesionDistanceCart
    if (distance_scale_factor > 1.0) {
      distance_scale_factor = 1.0;
    }

    // It tries to attach the CAR-T cell to the tumor cell with probability
    // kAdhesionRateCart * oncoprotein_scale_factor * distance_scale_factor *
    // kDtMechanics
    if (rng->Uniform(0.0, 1.0) < kAdhesionRateCart * oncoprotein_scale_factor *
                                     distance_scale_factor * kDtMechanics) {
// avoid race condition. Only one cell can be attached to the tumor cell.
#pragma omp critical
      {
        // We need to check again if the victim is not attached to a CAR-T cell
        // yet. This could be made more efficiently with a semaphore for each
        // cancer cell
        if (!victim->IsAttachedToCart()) {
          attached_to_tumor_cell_ = true;
          attached_cell_ptr_ = victim->GetAgentPtr<TumorCell>();
          victim->SetAttachedToCart(true);
        }
      }
    }
  }
}

// Try to induce apoptosis
bool CarTCell::TryToInduceApoptosis(bdm::AgentPointer<TumorCell> attached_cell,
                                    Random* rng) const {
  // If there is no attached cell (this should not happen)

  if (!attached_to_tumor_cell_) {
    return false;
  }

  // factor of how high is the oncoprotein levelof the cancer cell
  real_t scale_factor =
      (attached_cell->GetOncoproteinLevel() - kOncoproteinLimit) /
      kOncoproteinDifference;
  // Clamp scale_factor to be in [0,1]
  if (scale_factor > 1.0) {
    scale_factor = 1.0;
  }
  // If oncoprotein level is lower than the limit the cancer cell does not
  // become apoptotic
  if (scale_factor < 0.0) {
    // scale_factor = 0.0; the probability is going to be 0 so return the
    // function is the most efficient
    return false;
  }
  // CAR-T success of killing probability: aggressive cancer cells (high
  // oncoprotein level) are more likely to be killed
  const bool succeeded =
      rng->Uniform(0.0, 1.0) < kKillRateCart * scale_factor * kDtMechanics;

  // The CAR-T has succeeded to induce apoptosis on the Cancer Cell
  if (succeeded) {
    attached_cell->StartApoptosis();
  }

  return succeeded;
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
    // CAR-T do not change immunostimulatory factor levels
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
  Simulation* sim = Simulation::GetActive();
  // Run only every kDtCycle minutes, fmod does not work with the type
  // returned by GetSimulatedTime()
  if (sim->GetScheduler()->GetSimulatedSteps() % kStepsPerCycle != 0) {
    return;
  }

  if (auto* cell = dynamic_cast<CarTCell*>(agent)) {
    switch (cell->GetState()) {
      case CarTCellState::kAlive: {
        // the cell is growing to real_t its size before mitosis
        // Probability of death= 1/CurrentLiveTime,avoiding division by 0
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
            cell->GetAttachedCellPointer()->SetAttachedToCart(false);
            cell->SetAttachedCellPointer(nullptr);
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
          ExecutionContext* ctxt = sim->GetExecutionContext();
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
