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

#ifndef CORE_UTIL_UTILS_AUX_H_
#define CORE_UTIL_UTILS_AUX_H_

#include "core/container/math_array.h"
#include "core/operation/operation.h"
#include "core/operation/operation_registry.h"
#include "core/real_t.h"
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

namespace bdm {
/// Forward declaration of TumorCell class
class TumorCell;

/// Sample a positive Gaussian value
///
/// Samples a Gaussian value with given mean and standard deviation.
/// All negative values are mapped to zero to ensure positive results.
///
/// @param mean Mean value of the Gaussian distribution
/// @param sigma Standard deviation of the Gaussian distribution
/// @return Sampled positive value (negative values mapped to zero)
real_t SamplePositiveGaussian(float mean, float sigma);

/// Create a spherical arrangement of tumor cells
///
/// Generates a vector of 3D positions for tumor cells arranged in a spherical
/// pattern with the specified radius. The cells are positioned to form a
/// forces-stable structure. The achieved tumor mimics a heterogeneous organoid
/// thanks to the different cancer cell types that can be randomly generated
///
/// @param sphere_radius Radius of the spherical tumor in micrometers
/// @return Vector of 3D positions where tumor cells should be placed
std::vector<Real3> CreateSphereOfTumorCells(real_t sphere_radius);

/// Compute tumor statistics and characteristics
///
/// Analyzes the current tumor population to compute the number of tumor cells
/// of each type and the overall radius of the tumor mass. In addition, it
/// computes the average oncoprotein level and oxygen across all tumor cells.
///
/// @return Tuple containing:
///   - Total number of tumor cells
///   - Number of type 1 tumor cells (most aggressive)
///   - Number of type 2 tumor cells
///   - Number of type 3 tumor cells
///   - Number of type 4 tumor cells (least aggressive)
///   - Number of type 5 tumor cells (dead)
///   - Number of living CAR-T cells (not apoptotic)
///   - Current tumor radius in micrometers
///   - Average oncoprotein level across all tumor cells
///   - Average oxygen level across all tumor cells
std::tuple<size_t, size_t, size_t, size_t, size_t, size_t, size_t, real_t,
           real_t, real_t>
AnalyzeTumor();

/// Generates a random direction unitary vector
///
/// @return A 3D vector representing a random direction
Real3 GenerateRandomDirection();

/// Spawns a Dosage of CAR-T cells around the tumor
///
/// Called automatically by the simulation scheduler at the specified frequency.
/// Spawns a dosage of CAR-T cells around the tumor following the map
/// gKTreatment where
///   - The key represents the day of treatment (starting from day 0).
///   - The value represents the number of CAR-T cells administered on that day.
struct SpawnCart : public StandaloneOperationImpl {
  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  BDM_OP_HEADER(SpawnCart);

 public:
  void SetFrequency(uint64_t frequency) { frequency_ = frequency; }
  [[nodiscard]] uint64_t GetFrequency() const { return frequency_; }

 private:
  /// Frequency of output (every N simulation steps)
  uint64_t frequency_ = 1;

  void operator()() override;
};

/// Register SpawnCart operation with CPU as compute target
// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
inline BDM_REGISTER_OP(SpawnCart, "SpawnCart", kCpu);

/// Operation for outputting simulation summary data to a CSV file
///
/// This operation collects and outputs summary statistics about the simulation
/// state to output/final_data.csv for post-processing and analysis. It includes
/// information about cell populations, tumor characteristics, and other
/// relevant metrics.
struct OutputSummary : public StandaloneOperationImpl {
  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  BDM_OP_HEADER(OutputSummary);

 public:
  void SetFrequency(uint64_t frequency) { frequency_ = frequency; }
  [[nodiscard]] uint64_t GetFrequency() const { return frequency_; }

 private:
  /// Frequency of output (every N simulation steps)
  uint64_t frequency_ = 1;

  void operator()() override;
};

/// Register OutputSummary operation with CPU as compute target
// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
inline BDM_REGISTER_OP(OutputSummary, "OutputSummary", kCpu);

}  // namespace bdm

#endif  // CORE_UTIL_UTILS_AUX_H_
