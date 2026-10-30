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
#ifndef CART_TUMOR_H_
#define CART_TUMOR_H_

#include "params/hyperparams.h"
#include "utils/utils_aux.h"

namespace bdm {

/// List the diffused substances
enum Substances { kImmunostimulatoryFactor, kOxygen, kGlucose };

/// Run one simulation with the given parameters. It can be called several
/// times in the same process (e.g. from Python through cppyy)
/// @param input_params Parameters of the simulation
/// @param result If not null, it is filled in memory with the same data as
/// the output CSV files
/// @param write_csv Write the output CSV files in output/
/// @param verbose Print the parameters and the progress messages
/// @param allow_visualization Keep the visualization export of bdm.toml and
/// the command line, otherwise it is disabled for this run
/// @param argc Number of BioDynaMo command line arguments in argv
/// @param argv BioDynaMo command line arguments (nullptr for none)
/// @return 0 if the simulation finished, 1 if there was an error
int RunSimulation(const SimParam& input_params,
                  SimulationResult* result = nullptr, bool write_csv = false,
                  bool verbose = false, bool allow_visualization = false,
                  int argc = 0, const char** argv = nullptr);

/// Function declaration for the main simulation
int Simulate(int argc, const char** argv);

}  // namespace bdm

#endif  // CART_TUMOR_H_
