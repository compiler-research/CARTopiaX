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

#ifndef CORE_UTIL_UTILS_AUX_H_
#define CORE_UTIL_UTILS_AUX_H_

#include "biodynamo.h"
#include "core/util/random.h"
#include "hyperparams.h"
#include "tumor_cell.h"

namespace bdm {
    /// Forward declaration of TumorCell class
    class TumorCell;

/** @brief Sample a positive Gaussian value
 * 
 * Samples a Gaussian value with given mean and standard deviation.
 * All negative values are mapped to zero to ensure positive results.
 * 
 * @param mean Mean value of the Gaussian distribution
 * @param sigma Standard deviation of the Gaussian distribution
 * @return Sampled positive value (negative values mapped to zero)
 */
real_t SamplePositiveGaussian(float mean, float sigma);

/** @brief Create a spherical arrangement of tumor cells
 * 
 * Generates a vector of 3D positions for tumor cells arranged in a spherical
 * pattern with the specified radius. The cells are positioned to form a
 * realistic tumor structure.
 * 
 * @param sphere_radius Radius of the spherical tumor in micrometers
 * @return Vector of 3D positions where tumor cells should be placed
 */
std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius);

/** @brief Compute tumor statistics and characteristics
 * 
 * Analyzes the current tumor population to compute the number of tumor cells
 * of each type and the overall radius of the tumor mass.
 * 
 * @return Tuple containing:
 *   - Number of type 1 tumor cells (most aggressive)
 *   - Number of type 2 tumor cells
 *   - Number of type 3 tumor cells  
 *   - Number of type 4 tumor cells (least aggressive)
 *   - Number of type 5 tumor cells (dead)
 *   - Total number of tumor cells
 *   - Current tumor radius in micrometers
 */
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, real_t> ComputeNumberTumorCellsAndRadius();

/** @brief Operation for outputting simulation summary data to CSV files
 * 
 * This operation collects and outputs summary statistics about the simulation
 * state to CSV files for post-processing and analysis. It includes information
 * about cell populations, tumor characteristics, and other relevant metrics.
 */
struct OutputSummary : public StandaloneOperationImpl {
  BDM_OP_HEADER(OutputSummary);

  ///Frequency of output (every N simulation steps)
  uint64_t frequency_ = 1;

  ///Collects current simulation data and writes it to CSV files: Called automatically by the simulation scheduler at the specified frequency.
  void operator()() override;
};

/// Register OutputSummary operation with CPU as compute target
inline BDM_REGISTER_OP(OutputSummary, "OutputSummary", kCpu);


}  // namespace bdm

#endif  // CORE_UTIL_UTILS_AUX_H_