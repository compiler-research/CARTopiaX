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

namespace bdm {

/// List the diffused substances
enum Substances { kImmunostimulatoryFactor, kOxygen };

/// Function declaration for the main simulation
int Simulate(int argc, const char** argv);

}  // namespace bdm

#endif  // CART_TUMOR_H_