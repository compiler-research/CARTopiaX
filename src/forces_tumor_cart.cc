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

#include "forces_tumor_cart.h"


namespace bdm {

Real4 InteractionVelocity::Calculate(const Agent* lhs, const Agent* rhs) const {

  auto* a = dynamic_cast<const Cell*>(lhs);
  auto* b = dynamic_cast<const Cell*>(rhs);

  // Ignore self-interaction
  if (a->GetUid() == b->GetUid())
    return {0.0, 0.0, 0.0, 0.0};

  Real3 displacement = a->GetPosition() - b->GetPosition();

  //For periodic boundary conditions, we need to adjust the displacement
  displacement[0] = displacement[0] - (kBoundedSpaceLength)*round(displacement[0]/(kBoundedSpaceLength));
  displacement[1] = displacement[1] - (kBoundedSpaceLength)*round(displacement[1]/(kBoundedSpaceLength));
  displacement[2] = displacement[2] - (kBoundedSpaceLength)*round(displacement[2]/(kBoundedSpaceLength));

  double dist_sq = displacement[0] * displacement[0] +
                   displacement[1] * displacement[1] +
                   displacement[2] * displacement[2];
  double distance = std::max(std::sqrt(dist_sq), 1e-5);

  double radius_a = a->GetDiameter() / 2.0;
  double radius_b = b->GetDiameter() / 2.0;
  double R = radius_a + radius_b;
  // R=16.8254;//Debug
  // std::cout << "Debug: R = " << R << ", distance = " << distance << std::endl;// Debug output
  double temp_r = 0.0;

  const TumorCell* a_tumor = dynamic_cast<const TumorCell*>(a);
  const TumorCell* b_tumor = dynamic_cast<const TumorCell*>(b);

  if (distance < R) {

    // 1 - d/R
    temp_r = 1.0 - distance / R; 
    // (1 - d/R)^2
    temp_r *= temp_r;

    double repulsion;
  // std::cout << "temp_r = " << temp_r<< std::endl;// Debug output

    
    if (a_tumor && b_tumor) {// two tumor cells
      repulsion = kRepulsionTumorTumor;//std::sqrt(kRepulsionTumorTumor * kRepulsionTumorTumor);
    } else if (!a_tumor && !b_tumor) {// two CAR-T cells
      repulsion = kRepulsionCartCart;//std::sqrt(kRepulsionCartCart*kRepulsionCartCart);
    } else {// one tumor cell and one CAR-T
      repulsion = std::sqrt(kRepulsionCartTumor *
                            kRepulsionTumorCart);
    }
  // std::cout << "repulsion = " << repulsion<< std::endl;// Debug output

    temp_r *= repulsion;
  }

  // std::cout << "temp_r after repulsion = " << temp_r<< std::endl;// Debug output


  // Adhesion
  double max_interaction_distance = kMaxRelativeAdhesionDistance * R;
  // max_interaction_distance=21.0318;//Debug
  // std::cout << "max_interaction_distance = " << max_interaction_distance << std::endl;// Debug output


  if (distance < max_interaction_distance) {
    // 1 - d/S
    double temp_a = 1.0 - distance / max_interaction_distance; 
    // (1-d/S)^2
    temp_a *= temp_a;

    // std::cout << "temp_a = " << temp_a << std::endl;// Debug output


    double adhesion;
    if (a_tumor && b_tumor) {// two tumor cells
      adhesion = kAdhesionTumorTumor;
    } else if (!a_tumor && !b_tumor) {// two CAR-T cells
      adhesion = kAdhesionCartCart;
    } else {// one tumor cell and one CAR-T
      adhesion = std::sqrt(kAdhesionCartTumor *
                            kAdhesionTumorCart);
    }

    // std::cout << "adhesion = " << adhesion << std::endl;// Debug output


    temp_a *= adhesion;
    temp_r -= temp_a;

    // std::cout << "temp_a after adhesion= " << temp_a << std::endl;// Debug output

  }

  if (std::abs(temp_r) < 1e-16) {
    return {0.0, 0.0, 0.0, 0.0};
  }
  double force_magnitude = temp_r / distance;



      //Debug Output volcities
    // std::ofstream file("output/intercation_velocities.csv", std::ios::app);
    // if (file.is_open()) {
      
    //   double total_minutes = Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
    //   Real3 position = a->GetPosition();
    //   // Write data to CSV file
    //   file << total_minutes << ",position"
    //    << position[0] << ","
    //    << position[1] << ","
    //    << position[2] << ",displacement"
    //     << displacement[0] << ","
    //     << displacement[1] << ","
    //     << displacement[2] << ",distance"
    //     << distance << ",force_magnitude"
    //     << force_magnitude << ",temp_r"
    //    << temp_r << "\n";
    // }
    // End Debug Output


  // return{0.,0.,0.,0.};//debug
  

  return {2*force_magnitude * displacement[0],
          2*force_magnitude * displacement[1],
          2*force_magnitude * displacement[2],
          0.0};  // 4th component is unused
}

InteractionForce* InteractionVelocity::NewCopy() const {
  return new InteractionVelocity();
}

}  // namespace bdm