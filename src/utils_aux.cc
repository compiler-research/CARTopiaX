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

#include "utils_aux.h"

namespace bdm {

// Samples a Gaussian value with given mean and standard deviation but all negative values are mapped to zero
real_t SamplePositiveGaussian(float mean, float sigma) {
  auto* random = Simulation::GetActive()->GetRandom();
  real_t value = random->Gaus(mean, sigma);
  if(value < 0.) {value = 0.;}
  return value;
}

// In a future version, consider using this Behavior but controlling the time in which it is applied so that it is executed always after the diffusion module
// ConsumptionSecretion::ConsumptionSecretion(const std::string& substance, real_t quantity_consumption, real_t quantity_secretion, real_t substance_saturation_density){
//   dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid(substance);
//   substance_saturation_density_ = substance_saturation_density;
//   this->SetQuantities(quantity_consumption,quantity_secretion);
//   AlwaysCopyToNew();
// }

// ConsumptionSecretion::ConsumptionSecretion(DiffusionGrid* dgrid, real_t quantity_consumption, real_t quantity_secretion,real_t substance_saturation_density){
//       dgrid_ = dgrid;
//       substance_saturation_density_ = substance_saturation_density;
//       this->SetQuantities(quantity_consumption,quantity_secretion);
//       AlwaysCopyToNew(); 
//     }

// void ConsumptionSecretion::Initialize(const NewAgentEvent& event) {
//   Base::Initialize(event);
//   auto* other = bdm_static_cast<ConsumptionSecretion*>(event.existing_behavior);
//   dgrid_ = other->dgrid_;
//   quantity_consumption_ = other->quantity_consumption_;
//   quantity_secretion_ = other->quantity_secretion_;
//   substance_saturation_density_ = other->substance_saturation_density_;
//   constant1_ = other->constant1_;
//   constant2_ = other->constant2_;
// }

// void ConsumptionSecretion::Run(Agent* agent) {

//   const auto& pos = agent->GetPosition();

//   // // std::cout<<"constant1: "<<constant1_<<" constant2: "<<constant2_<<std::endl; //Debug
//   // //Debug
//   // constant1_=0.;
//   // constant2_=1.3;
//   real_t conc;
//   real_t new_conc;
//   #pragma omp critical
//   {
//     conc = dgrid_->GetValue(pos);
//     new_conc = (conc + constant1_) / constant2_; // Equation solution
//     dgrid_->ChangeConcentrationBy(pos, new_conc - conc, InteractionMode::kAdditive, false);//CHANGE to use SetConcentration from DiffusionThomasAlgorithm
//   }
    
//     // double current_time = Simulation::GetActive()->GetScheduler()->GetSimulatedSteps()* kDt; // Get the current time step in minutes
//     // std::ofstream file("output/consumptions_mine.csv", std::ios::app);////Debug
//     // if (file.is_open()) {
//     // file  << current_time << "," << conc << "," << constant1_ << "," << constant2_ << "," << new_conc <<"\n";
//     // }
//     // double current_time = Simulation::GetActive()->GetScheduler()->GetSimulatedSteps()* kDt; // Get the current time step in minutes
//     // std::ofstream file("output/consumptions_mine" + std::to_string(current_time) + ".csv", std::ios::app);////Debug//(12*60)) + ".csv", std::ios::app);
//     // if (file.is_open()) {
//     // file  << conc << "," << constant1_ << "," << constant2_ << "," << new_conc <<"\n";


//   // }
//   // Debug: print all parameters and values if new concentration > 100
//   // if (new_conc > 100) {
//   //   std::cout << "Debug ConsumptionSecretion:\n";
//   //   std::cout << "  Position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
//   //   std::cout << "  conc (current concentration): " << conc << "\n";
//   //   std::cout << "  constant1: " << constant1_ << "\n";
//   //   std::cout << "  constant2: " << constant2_ << "\n";
//   //   std::cout << "  new_conc: " << new_conc << "\n";
//   //   std::cout << "  quantity_consumption_: " << quantity_consumption_ << "\n";
//   //   std::cout << "  quantity_secretion_: " << quantity_secretion_ << "\n";
//   //   std::cout << "  substance_saturation_density_: " << substance_saturation_density_ << "\n";
//   //   std::cout << "  kDt: " << kDt << "\n";
//   //   std::cout << "  kDefaultVolumeNewTumorCell: " << kDefaultVolumeNewTumorCell << "\n";
//   //   std::cout << "  kVoxelVolume: " << kVoxelVolume << "\n";
//   //   //stop simulation
//   //   throw std::runtime_error("Aborting simulation due to high concentration in ConsumptionSecretion.");
//   // }
// }

// void ConsumptionSecretion::SetQuantities(real_t quantity_consumption, real_t quantity_secretion) {// Set quantities for consumption and secretion by giving them already scaled
//   quantity_consumption_ = quantity_consumption;
//   quantity_secretion_ = quantity_secretion;
// }

// void ConsumptionSecretion::ComputeConstants(real_t total_volume) {
//   //compute the constants for the differential equation explicit solution
//   //dt*(cell_volume/voxel_volume)*quantity_secretion*substance_saturation =  dt · (V_k / V_voxel) · S_k · ρ*_k) 
//   constant1_ = quantity_secretion_ * substance_saturation_density_ * kDt * (total_volume / kVoxelVolume);// Scale by the volume of the cell in the Voxel and time step
//   //1 + dt*(cell_volume/voxel_volume)*(quantity_secretion + quantity_consumption ) = [1 + dt · (V_k / V_voxel) · (S_k + U_k)]
//   constant2_ = 1 + kDt * (total_volume/ kVoxelVolume) * (quantity_secretion_ + quantity_consumption_);// Scale by the volume of the cell in the Voxel and time step
// }

std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius) {
    // V = (4/3)*pi*r^3 = (pi/6)*diameter^3
    real_t cell_radius = std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi)/2;

    std::vector<Real3> positions;

    // Hexagonal close-packing spacing
    real_t spacing_x = cell_radius * std::sqrt(3.0);
    real_t spacing_y = cell_radius * 2.0;
    real_t spacing_z = cell_radius * std::sqrt(3.0);

    int zc = 0;
    for (real_t z = -sphere_radius; z < sphere_radius; z += spacing_z, ++zc) {
        int xc = 0;
        for (real_t x = -sphere_radius; x < sphere_radius; x += spacing_x, ++xc) {
            int yc = 0;
            for (real_t y = -sphere_radius; y < sphere_radius; y += spacing_y, ++yc) {

                // Compute cell center with HCP offset
                real_t px = x + (zc % 2) * 0.5 * cell_radius;
                real_t py = y + (xc % 2) * cell_radius;
                real_t pz = z;

                real_t dist = std::sqrt(px * px + py * py + pz * pz);

                if (dist <= sphere_radius) {
                    positions.push_back({px, py, pz});
                }
            }
        }
    }

    return positions;
}

//Function to compute the number of tumor cells of each type and the radius of the tumor
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, real_t> ComputeNumberTumorCellsAndRadius() {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  size_t total_num_tumor_cells = 0;
  size_t num_tumor_cells_type1 = 0;
  size_t num_tumor_cells_type2 = 0;
  size_t num_tumor_cells_type3 = 0;
  size_t num_tumor_cells_type4 = 0;
  size_t num_tumor_cells_type5_dead = 0;

  real_t max_dist_sq = 0.0;

  rm->ForEachAgent([&](const Agent* agent) {
    if (auto* tumor_cell = dynamic_cast<const TumorCell*>(agent)) {
      total_num_tumor_cells++;
      const auto& pos = agent->GetPosition();
      real_t dist_sq = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
      if (dist_sq > max_dist_sq) {
        max_dist_sq = dist_sq;
      }

      // Count tumor cells by type
      switch (tumor_cell->GetType()) {
        case 1: num_tumor_cells_type1++; break;
        case 2: num_tumor_cells_type2++; break;
        case 3: num_tumor_cells_type3++; break;
        case 4: num_tumor_cells_type4++; break;
        case 5: num_tumor_cells_type5_dead++; break;
        default: break;
      }
    }
  });
  return {total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, std::sqrt(max_dist_sq)};
}

// Function to output summary CSV
void OutputSummary::operator()() {
  auto* scheduler = Simulation::GetActive()->GetScheduler();
  auto current_step = scheduler->GetSimulatedSteps();

  if (current_step % frequency_ == 0) {
    std::ofstream file("output/final_data.csv", std::ios::app);
    if (file.is_open()) {
      if (current_step == 0) {
        file << "total_days,total_hours,total_minutes,tumor_radius,num_cells,num_tumor_cells,tumor_cells_type1,tumor_cells_type2,tumor_cells_type3,tumor_cells_type4,tumor_cells_type5_dead\n";// Header for CSV file
      }

      // Calculate time in days, hours, minutes
      double total_minutes = Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
      double total_hours = total_minutes / 60.0;
      double total_days = total_hours / 24.0;

      // Count total cells, tumor cells of each type and tumor radius
      size_t total_num_tumor_cells;
      size_t num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead;
      real_t tumor_radius;
      std::tie(total_num_tumor_cells, num_tumor_cells_type1, num_tumor_cells_type2, num_tumor_cells_type3, num_tumor_cells_type4, num_tumor_cells_type5_dead, tumor_radius) = ComputeNumberTumorCellsAndRadius();
      // Write data to CSV file
      file << total_days << ","
       << total_hours << ","
       << total_minutes << ","
       << tumor_radius << ","
       << Simulation::GetActive()->GetResourceManager()->GetNumAgents() << ","//number of cells
       << total_num_tumor_cells << ","
       << num_tumor_cells_type1 << ","
       << num_tumor_cells_type2 << ","
       << num_tumor_cells_type3 << ","
       << num_tumor_cells_type4 << ","
       << num_tumor_cells_type5_dead << "\n";
    }
  }
}

}  // namespace bdm