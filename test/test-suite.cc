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

#include "biodynamo.h"
#include <gtest/gtest.h>

// Googletest in combination with the provided CMakeLists.txt allows you to
// define tests in arbitrary .cc files in the `test/` folder. This file should
// serve as an inspiration for testing user-defined, custom behaviors, basic as
// well as compicated functions, or similar things. If you wish to add tests for
// specific aspects, you can either add them to the existing test-suite.cc file
// or create a new *.cc file in the `test/` folder. CMake will handle it
// automatically. For more information regarding testing with Googletest, please
// consider the following sources:
// * https://google.github.io/googletest/primer.html
// * https://github.com/google/googletest

#define TEST_NAME typeid(*this).name()

namespace bdm {

// A function to test
int Compute42() { return 6 * 7; };

// Show how to compare two numbers
TEST(UtilTest, NumberTest) {
  // Expect equality
  EXPECT_EQ(Compute42(), 42);
}

// Test if we can add agents to the simulation
TEST(AgentTest, AddAgentsToSimulation) {
  // Create simulation
  Simulation simulation(TEST_NAME);

  // Add some cells to the simulation
  auto* rm = simulation.GetResourceManager();
  uint8_t expected_no_cells{20};
  for (int i = 0; i < expected_no_cells; i++) {
    auto* cell = new Cell(30);
    rm->AddAgent(cell);
  }

  // Test if all 20 cells are in the simulation
  auto no_cells = rm->GetNumAgents();
  EXPECT_EQ(expected_no_cells, no_cells);
}

}  // namespace bdm
