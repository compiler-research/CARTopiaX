/*
 * Copyright 2026 compiler-research.org
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
 * This file contains a model developed for the compiler-research.org 
 * organization.
 */

// Tests of the in-process entry points used by the Python driver
// (abm_calibration/optimize.py): SimParam::LoadParamsFromString and
// RunSimulation with the output stored in memory

#include "params/hyperparams.h"
#include "utils/utils_aux.h"
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <string>
#include <system_error>
#include <utility>
#include <vector>
#include "cart_tumor.h"

namespace bdm {

namespace {

// Parameters with values different from the defaults, including some that
// change the values derived from them (dt_mechanics, tumor_shape)
constexpr const char* kCustomParams = R"({
  "seed": 7,
  "total_minutes_to_simulate": 600,
  "tumor_shape": "cylinder",
  "bound_space_toplogy": "closed",
  "initial_oxygen_level": 33.5,
  "dt_mechanics": 0.5,
  "treatment": {"2": 10},
  "output_information_dependent_on_radius": true,
  "num_radius_intervals": 4
})";

// Short simulation with a small tumor to keep the tests fast: 10 minutes with
// dt_step = 0.1 (100 steps) and an output every 50 steps (minutes 0, 5, 10)
constexpr const char* kShortSimulation = R"({
  "seed": 7,
  "total_minutes_to_simulate": 10,
  "initial_spherical_tumor_radius": 30,
  "resolution_grid_substances": 20,
  "output_csv_interval": 50,
  "treatment": {"0": 5}
})";

// Number of outputs of kShortSimulation
constexpr size_t kShortSimulationOutputs = 3;

constexpr const char* kFinalDataCsv = "output/final_data.csv";

// Expect the same loaded and derived parameters
void ExpectSameParams(const SimParam& a, const SimParam& b) {
  EXPECT_EQ(a.seed, b.seed);
  EXPECT_EQ(a.total_minutes_to_simulate, b.total_minutes_to_simulate);
  EXPECT_EQ(a.tumor_shape, b.tumor_shape);
  EXPECT_EQ(a.bound_space_toplogy, b.bound_space_toplogy);
  EXPECT_EQ(a.initial_oxygen_level, b.initial_oxygen_level);
  EXPECT_EQ(a.treatment, b.treatment);
  EXPECT_EQ(a.output_information_dependent_on_radius,
            b.output_information_dependent_on_radius);
  EXPECT_EQ(a.num_radius_intervals, b.num_radius_intervals);
  EXPECT_EQ(a.dt_step, b.dt_step);
  EXPECT_EQ(a.output_csv_interval, b.output_csv_interval);
  EXPECT_EQ(a.steps_in_one_day, b.steps_in_one_day);
  EXPECT_EQ(a.steps_per_cell_cycle, b.steps_per_cell_cycle);
  EXPECT_EQ(a.diffuse_glucose_on_z_axis, b.diffuse_glucose_on_z_axis);
  EXPECT_EQ(a.voxel_volume, b.voxel_volume);
  EXPECT_EQ(a.oncoprotein_difference, b.oncoprotein_difference);
  EXPECT_EQ(a.max_squared_distance_cart_moving_towards_tumor_cell,
            b.max_squared_distance_cart_moving_towards_tumor_cell);
}

// Whether a file exists and its last modification time, to check that a
// simulation does not write it
std::pair<bool, std::filesystem::file_time_type> GetFileState(
    const std::string& path) {
  std::error_code error;
  const auto time = std::filesystem::last_write_time(path, error);
  if (error) {
    return {false, {}};
  }
  return {true, time};
}

}  // namespace

// case: when a JSON string loads the same parameters as a file with that
// content
TEST(LoadParamsFromStringTest, SameAsLoadParamsFromFile) {
  const auto path =
      std::filesystem::temp_directory_path() / "cartopiax_test_params.json";
  {
    std::ofstream file(path);
    file << kCustomParams;
  }
  SimParam from_file;
  from_file.LoadParams(path.string());
  std::filesystem::remove(path);

  SimParam from_string;
  from_string.LoadParamsFromString(kCustomParams);

  ExpectSameParams(from_file, from_string);
}

// case: the values in the string are loaded and the derived ones recomputed
TEST(LoadParamsFromStringTest, LoadsValuesAndDerivedParams) {
  SimParam params;
  params.LoadParamsFromString(kCustomParams);

  EXPECT_EQ(params.seed, 7);
  EXPECT_EQ(params.total_minutes_to_simulate, 600);
  EXPECT_EQ(params.tumor_shape, TumorShape::kCylinder);
  EXPECT_EQ(params.bound_space_toplogy, Param::BoundSpaceMode::kClosed);
  EXPECT_DOUBLE_EQ(params.initial_oxygen_level, 33.5);
  // The treatment replaces the default one
  EXPECT_EQ(params.treatment, (std::map<int, int>{{2, 10}}));
  EXPECT_TRUE(params.output_information_dependent_on_radius);
  EXPECT_EQ(params.num_radius_intervals, 4);

  // Derived from dt_mechanics = 0.5 and tumor_shape = cylinder
  EXPECT_DOUBLE_EQ(params.dt_step, 0.5);
  EXPECT_EQ(params.output_csv_interval, 12 * 60 * 2);
  EXPECT_EQ(params.steps_in_one_day, static_cast<size_t>(24 * 60 * 2));
  EXPECT_EQ(params.steps_per_cell_cycle, 12);
  EXPECT_FALSE(params.diffuse_glucose_on_z_axis);
}

// case: an empty or invalid string keeps the default parameters
TEST(LoadParamsFromStringTest, EmptyOrInvalidStringKeepsDefaults) {
  SimParam defaults;
  defaults.LoadParamsFromString("{}");

  SimParam from_empty;
  from_empty.LoadParamsFromString("");
  ExpectSameParams(defaults, from_empty);

  SimParam from_invalid;
  from_invalid.LoadParamsFromString("{ not json");
  ExpectSameParams(defaults, from_invalid);
}

// case: the output is stored in memory with the columns of final_data.csv and
// the CSV file is not written
TEST(RunSimulationTest, FillsFinalDataInMemory) {
  const auto csv_state = GetFileState(kFinalDataCsv);

  SimParam params;
  params.LoadParamsFromString(kShortSimulation);
  SimulationResult result;
  ASSERT_EQ(RunSimulation(params, &result), 0);

  // Same columns as the header of output/final_data.csv
  const std::vector<std::string> columns = {
      "total_days",
      "total_hours",
      "total_minutes",
      "tumor_radius",
      "num_cells",
      "num_tumor_cells",
      "tumor_cells_type1",
      "tumor_cells_type2",
      "tumor_cells_type3",
      "tumor_cells_type4",
      "tumor_cells_type5_dead",
      "dead_tumor_cells_lack_of_oxygen",
      "dead_tumor_cells_lack_of_glucose",
      "dead_tumor_cells_random_natural_causes",
      "dead_tumor_cells_cartcell_kill",
      "num_alive_cart",
      "num_dead_cart",
      "average_oncoprotein",
      "average_oxygen_cancer_cells",
      "average_oxygen_all_cells",
      "average_glucose_cancer_cells",
      "average_glucose_all_cells",
      "average_radius_distance_living_cart_cells"};
  ASSERT_EQ(result.final_data.size(), columns.size());
  for (const auto& column : columns) {
    ASSERT_EQ(result.final_data.count(column), 1U) << column;
    EXPECT_EQ(result.final_data.at(column).size(), kShortSimulationOutputs)
        << column;
  }

  const auto& minutes = result.final_data.at("total_minutes");
  EXPECT_DOUBLE_EQ(minutes[0], 0);
  EXPECT_DOUBLE_EQ(minutes[1], 5);
  EXPECT_DOUBLE_EQ(minutes[2], 10);
  EXPECT_GT(result.final_data.at("num_tumor_cells")[0], 0);
  // The CAR-T cells of the dose of day 0 are counted at minute 0
  EXPECT_DOUBLE_EQ(result.final_data.at("num_alive_cart")[0], 5);

  // Without output_information_dependent_on_radius there is no radial data
  EXPECT_TRUE(result.radius_data.empty());
  // The CSV file is not written
  EXPECT_EQ(GetFileState(kFinalDataCsv), csv_state);
}

// case: the radial output is stored in memory with the CSV column names
TEST(RunSimulationTest, FillsRadiusDataInMemory) {
  auto config = nlohmann::json::parse(kShortSimulation);
  config["output_information_dependent_on_radius"] = true;
  config["max_radius_analysis_csv_dependent_on_radius"] = 300;
  config["num_radius_intervals"] = 2;

  SimParam params;
  params.LoadParamsFromString(config.dump());
  SimulationResult result;
  ASSERT_EQ(RunSimulation(params, &result), 0);

  // Time columns and 22 metrics for each of the 2 radius intervals
  EXPECT_EQ(result.radius_data.size(), 3U + 22U * 2U);
  for (const auto& [column, values] : result.radius_data) {
    EXPECT_EQ(values.size(), kShortSimulationOutputs) << column;
  }
  EXPECT_EQ(result.radius_data.count("tumor_cells_type5_dead_radius_0_to_150"),
            1U);
  EXPECT_EQ(result.radius_data.count("num_alive_cells_radius_150_to_300"), 1U);
}

// case: repeated runs in the same process with the same seed give the same
// output, so no state is kept between runs
TEST(RunSimulationTest, RepeatedRunsWithSameSeedGiveSameResult) {
  // With several threads the order in which agents are updated can change
  const int threads = omp_get_max_threads();
  omp_set_num_threads(1);

  SimParam params;
  params.LoadParamsFromString(kShortSimulation);
  SimulationResult first;
  SimulationResult second;
  EXPECT_EQ(RunSimulation(params, &first), 0);
  EXPECT_EQ(RunSimulation(params, &second), 0);

  omp_set_num_threads(threads);

  EXPECT_EQ(first.final_data, second.final_data);
}

// case: an unknown tumor shape returns an error code
TEST(RunSimulationTest, UnknownTumorShapeReturnsError) {
  SimParam params;
  params.LoadParamsFromString(kShortSimulation);
  params.tumor_shape = static_cast<TumorShape>(-1);
  SimulationResult result;
  EXPECT_EQ(RunSimulation(params, &result), 1);
  EXPECT_TRUE(result.final_data.empty());
}

}  // namespace bdm
