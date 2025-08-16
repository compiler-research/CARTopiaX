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

#include "cart_cell.h"

namespace bdm {

CartCell::CartCell(const Real3& position) {
  SetPosition(position);
  // Default state for new cells
  state_ = CartCellState::kAlive;
  // Initial timer_state for apoptotic state
  timer_state_ = 0;

  //volumes
  // Set default volume
  SetVolume(kDefaultVolumeNewCartCell);
  // Set default fluid fraction
  SetFluidFraction(kDefaultFractionFluidCartCell);
  // Set default nuclear volume
  SetNuclearVolume(kDefaultVolumeNucleusCartCell);


  // Pointer to oxygen diffusion grid
  oxygen_dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid("oxygen");
  // Pointer to immunostimulatory_factor diffusion grid
  immunostimulatory_factor_dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid("immunostimulatory_factor");
  // Initially not attached to a tumor cell
  attached_to_tumor_cell_ = false; 
  // Initialize attached cell pointer to null
  attached_cell_ = nullptr;

  // Initialize the velocity of the cell in the previous step to zero
  older_velocity_ = {0, 0, 0};


  SetCurrentLiveTime(kAverageMaximumTimeUntillApoptosisCart);

  //Add Consumption and Secretion
  // Set default oxygen consumption rate
  SetOxygenConsumptionRate(kDefaultOxygenConsumption);
  // Compute constants for all ConsumptionSecretion of Oxygen
  ComputeConstantsConsumptionSecretion();

}

// Cart cells can move if they are alive and not attached to a tumor cell
bool CartCell::DoesCellMove() {
  return (state_ == CartCellState::kAlive && !attached_to_tumor_cell_); 
}


real_t CartCell::GetTargetTotalVolume() {
  return GetTargetNucleusSolid() * (1 + GetTargetRelationCytoplasmNucleus()) / (1 - GetTargetFractionFluid());
}

// This method explicitly solves the system of exponential relaxation differential equation using a discrete 
// update step. It is used to shrink the volume (and proportions) smoothly toward a desired target
// volume over time whe the cell is apoptotic. The relaxations rate controls the speed of convergence
void CartCell::ChangeVolumeExponentialRelaxationEquation(real_t relaxation_rate_cytoplasm, real_t relaxation_rate_nucleus, real_t relaxation_rate_fluid) {
  // Exponential relaxation towards the target volume
  real_t current_total_volume = GetVolume();
  real_t fluid_fraction= GetFluidFraction();
  real_t nuclear_volume = GetNuclearVolume();

  real_t current_nuclear_solid = nuclear_volume * (1 - fluid_fraction);
  real_t current_cytoplasm_solid = (current_total_volume - nuclear_volume) * (1-fluid_fraction);

  real_t current_fluid = fluid_fraction * current_total_volume;

  // Update fluid volume
  real_t new_fluid = current_fluid + kDtCycle* relaxation_rate_fluid * (GetTargetFractionFluid() * current_total_volume - current_fluid);
  // Clamp to zero to prevent negative volumes
  if (new_fluid < 0.0) { new_fluid = 0.0; }
  
  real_t nuclear_fluid = new_fluid* ( nuclear_volume/ current_total_volume);
  // real_t cytoplasm_fluid = new_fluid - nuclear_fluid;

  real_t nuclear_solid = current_nuclear_solid + kDtCycle * relaxation_rate_nucleus * (GetTargetNucleusSolid() - current_nuclear_solid);
  // Clamp to zero to prevent negative volumes
  if (nuclear_solid < 0.0) { nuclear_solid = 0.0; }

  real_t target_cytoplasm_solid = GetTargetRelationCytoplasmNucleus() * GetTargetNucleusSolid();
  real_t cytoplasm_solid = current_cytoplasm_solid + kDtCycle * relaxation_rate_cytoplasm * (target_cytoplasm_solid - current_cytoplasm_solid);
  // Clamp to zero to prevent negative volumes
  if (cytoplasm_solid < 0.0) { cytoplasm_solid = 0.0; }

  real_t new_total_solid= nuclear_solid + cytoplasm_solid;

  real_t total_nuclear= nuclear_solid + nuclear_fluid;

  // real_t total_cytoplasm= cytoplasm_solid + cytoplasm_fluid;

  real_t new_volume = new_total_solid + new_fluid;

  // Avoid division by zero
  real_t new_fraction_fluid = new_fluid / (1e-16 + new_volume);
  
  // Update the cell's properties
  // if the volume has changed
  if (new_volume!= current_total_volume){
    SetVolume(new_volume);
    // Update constants for all ConsumptionSecretion of Oxygen and Immunostimulatory Factors
    ComputeConstantsConsumptionSecretion();
  }
  SetFluidFraction(new_fraction_fluid);
  SetNuclearVolume(total_nuclear);
}

//compute Displacement
Real3 CartCell::CalculateDisplacement(const InteractionForce* force,
                            real_t squared_radius, real_t dt) {

  // real_t h = dt;
  Real3 movement_at_next_step{0, 0, 0};
  // this should be chaged in a future version of BioDynaMo in order to have a cleaner code instead of hardcoding it here
  squared_radius=kSquaredMaxDistanceNeighborsForce;

  // the physics force to move the point mass
  Real3 translation_velocity_on_point_mass{0, 0, 0};

  // We check for every neighbor object if they touch us, i.e. push us
  // away and agreagate the velocities

  uint64_t non_zero_neighbor_forces = 0;
  if (!IsStatic()) {
    auto* ctxt = Simulation::GetActive()->GetExecutionContext();
    auto calculate_neighbor_forces =
        L2F([&](Agent* neighbor, real_t squared_distance) {
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
  // position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 * velocity(t - dt) ]
  movement_at_next_step += translation_velocity_on_point_mass * kDnew + older_velocity_ * kDold;

  //Debug
    //     std::ofstream file1("output/movement_at_next_step.csv", std::ios::app);
    // if (file1.is_open()) {

    //   // Calculate time in days, hours, minutes
    //   double total_minutes = Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
    //   // Write data to CSV file
    //   file1
    //    << "minute"<<total_minutes << ","
    //    <<"translation_velocity_on_point_mass"<< translation_velocity_on_point_mass[0]<< ","
    //    << translation_velocity_on_point_mass[1] << ","
    //    << translation_velocity_on_point_mass[2] << ","
    //    << "kDnew" << kDnew << ","
    //    << "older_velocity_" << older_velocity_[0] << ","
    //     << older_velocity_[1] << ","
    //     << older_velocity_[2] << ","
    //    << "kDold" << kDold << "\n";
    // }

  //Debug Output forces
    // std::ofstream file("output/forces.csv", std::ios::app);
    // if (file.is_open()) {
      
    //   // Calculate time in days, hours, minutes
    //   double total_minutes = Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
    //   double modulus_total_displacement = movement_at_next_step[0] * movement_at_next_step[0] +
    //                                      movement_at_next_step[1] * movement_at_next_step[1] +
    //                                      movement_at_next_step[2] * movement_at_next_step[2];
    //   modulus_total_displacement = std::sqrt(modulus_total_displacement);
    //   Real3 position = GetPosition();
    //   // Write data to CSV file
    //   file 
    //    << total_minutes << ","
    //     << position[0] << ","
    //    << position[1] << ","
    //    << position[2] << ","
    //    << movement_at_next_step[0] << ","
    //    << movement_at_next_step[1] << ","
    //    << movement_at_next_step[2] << ","
    //    << modulus_total_displacement << "\n";
    // }
    // End Debug Output

  older_velocity_ = translation_velocity_on_point_mass;

  // Displacement
  return movement_at_next_step;
}

//Compute new oxygen or immunostimulatory factor concentration after consumption/ secretion
real_t CartCell::ConsumeSecreteSubstance(int substance_id, real_t old_concentration) {
  real_t res;
  if (substance_id == oxygen_dgrid_->GetContinuumId()) {
    // consuming oxygen
    res= (old_concentration + constant1_oxygen_) / constant2_oxygen_;
  } else if (substance_id == immunostimulatory_factor_dgrid_->GetContinuumId()) {
    // This point should never be reached
    res= old_concentration;
  } else {
    throw std::invalid_argument("Unknown substance id: " + std::to_string(substance_id));
  }
  return res;
}

//Recompute Consumption constants whenever oxygen_consumption_rate_ or the volume changes
void CartCell::ComputeConstantsConsumptionSecretion() {

  // constant1_= dt · (V_k / V_voxel) · S_k · ρ*_k) 
  // constant2_ = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // where:
  // S_k    = secretion rate of cell k
  // U_k    = uptake (consumption) rate of cell k
  // ρ*_k   = saturation (target) density for secretion
  // V_k    = volume of the cell k
  // V_voxel = volume of the voxel containing the cell
  // dt     = simulation time step
  real_t volume = GetVolume();
  //compute the constants for the differential equation explicit solution: for oxygen and immunostimulatory factor
  //dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt · (V_k / V_voxel) · S_k · ρ*_k) 
  constant1_oxygen_ = 0.;
  //1 + dt*(cell_volume/voxel_volume)*(quantity_secretion + quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
  // Scale by the volume of the cell in the Voxel and time step
  constant2_oxygen_ = 1 + kDtSubstances * (volume/kVoxelVolume) * (oxygen_consumption_rate_);
}

/// Main behavior executed at each simulation step
void StateControlCart::Run(Agent* agent) {

  auto* sim = Simulation::GetActive();
  if(sim->GetScheduler()->GetSimulatedSteps() % kStepsPerCycle != 0){return;}// Run only every kDtCycle minutes, fmod does not work with the type returned by GetSimulatedTime()

  if (auto* cell = dynamic_cast<CartCell*>(agent)) {

    switch (cell->GetState())
    {
      case CartCellState::kAlive:{//the cell is growing to real_t its size before mitosis

        if (sim->GetRandom()->Uniform(1.0) < kDtCycle/std::max(cell->GetCurrentLiveTime(), 1e-10)) { // Probability of death= 1/CurrentLiveTime, avoiding division by 0
          //the cell Dies
          cell->SetState(CartCellState::kApoptotic);
          cell->SetTimerState(0);  // Reset timer_state, it should be 0 anyway
          // Set target volume to 0 (the cell will shrink)
          cell->SetTargetCytoplasmSolid(0.0); 
          cell->SetTargetNucleusSolid(0.0); 
          cell->SetTargetFractionFluid(0.0); 
          cell->SetTargetRelationCytoplasmNucleus(0.0);
          //Reduce oxygen consumption
          cell->SetOxygenConsumptionRate(cell->GetOxygenConsumptionRate()*kReductionConsumptionDeadCells);
          cell->ComputeConstantsConsumptionSecretion(); // Update constants for all Consumption of oxygen
          if (cell->IsAttachedToTumorCell()) {// Detach from tumor cell if it was attached
            cell->GetAttachedCell()->SetAttachedToCart(false);
            cell->SetAttachedCell(nullptr);
            cell->SetAttachedToTumorCell(false);
          }
        } else{
          // decrease current life time
          cell->SetCurrentLiveTime((cell->GetCurrentLiveTime() - (kDtCycle*kDtCycle)));
        }
        break;
      }
      case CartCellState::kApoptotic:{
        cell->SetTimerState(cell->GetTimerState() + kDtCycle); 

        cell->ChangeVolumeExponentialRelaxationEquation(kVolumeRelaxarionRateCytoplasmApoptotic,
                                                kVolumeRelaxarionRateNucleusApoptotic,
                                                kVolumeRelaxarionRateFluidApoptotic);
                                                
        if (kTimeApoptosis < cell->GetTimerState()) { // If the timer_state exceeds the time to transition (this is a fixed duration transition)
          //remove the cell from the simulation
          auto* ctxt = sim->GetExecutionContext();
          ctxt->RemoveAgent(agent->GetUid());
        }
        break;
      }
      default:{
        Log::Error("StateControlCart::Run", "Unknown CartCellState");
        break;
      }
    }
  } else {
    Log::Error("StateControlCart::Run", "SimObject is not a CartCell");
  }
}

}  // namespace bdm