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

#ifndef TUMOR_HYPERPARAMS_H_
#define TUMOR_HYPERPARAMS_H_

#include <map>
#include <string>

namespace bdm {



/** @name Tumor Cell Hyperparameters
 *  @brief Parameters controlling tumor cell behavior and properties. This should be in a not recompiled parameters file in the future
 *  @{
 */

/** @brief Rate of immunostimulatory factor secretion by tumor cells (per minute) */
constexpr real_t kRateSecretionImmunostimulatoryFactor= 10.0;

/** @brief Saturation density of immunostimulatory factor for tumor cells */
constexpr real_t kSaturationDensityImmunostimulatoryFactor = 1.0;

/** @brief Mean level of oncoprotein expression in tumor cells */
constexpr real_t kOncoproteinMean = 1.0;

/** @brief Standard deviation of oncoprotein expression in tumor cells */
constexpr real_t kOncoproteinStandardDeviation = 0.25;

/** @brief Oxygen saturation level required for tumor cell proliferation (mmHg) */
constexpr real_t kOxygenSaturationInProliferation = 38.0;

/** @brief Minimum oxygen level threshold for tumor cell proliferation (mmHg) */
constexpr real_t kOxygenLimitForProliferation = 10.0;

/** @brief Oxygen level threshold below which necrosis begins (mmHg) */
constexpr real_t kOxygenLimitForNecrosis = 5.0;

/** @brief Oxygen level for maximum necrosis probability (mmHg) */
constexpr real_t kOxygenLimitForNecrosisMaximum= 2.5;
/** @brief Transition rates for Ki67 expression states (commented out)
 * 
 * Transition rates for Ki67 expression states in min⁻¹:
 * - Rate from 0 to 1 (ki67 negative to ki67 positive pre-mitotic)
 * - Rate from 1 to 2 (ki67 positive pre-mitotic to ki67 positive post-mitotic)  
 * - Rate from 2 to 0 (ki67 positive post-mitotic to ki67 negative)
 */
// constexpr real_t kTransitionRateKi67[] = {// Transition rates for Ki67 expression states in min**-1
//     1.0 / (3.62 * 60.0),  // Rate from 0 to 1 (ki67 negative to ki67 positive pre-mitotic)
//     1.0 / (13.0 * 60.0),  // Rate from 1 to 2 (ki67 positive pre-mitotic to ki67 positive post-mitotic)
//     1.0 / (2.5 * 60.0)    // Rate from 2 to 0 (ki67 positive post-mitotic to ki67 negative)
// };

/** @brief Time until a lysed necrotic cell is removed from simulation (minutes) */
constexpr real_t kTimeLysis = 60*24*60.;

/** @brief Rate of cell division (min⁻¹) */
constexpr real_t kDivisionRate = 0.02717 / 60.0;

/** @brief Maximum necrosis rate for tumor cells in hypoxic conditions (min⁻¹) */
constexpr real_t kMaximumNecrosisRate= 1.0 / (6.0 * 60.0);

/** @brief Default oxygen consumption rate of tumor cells */
constexpr real_t kDefaultOxygenConsumption = 10.0;
/** @name Volume Parameters
 *  @brief Default volume parameters for tumor cells
 *  @{
 */

/** @brief Default total volume of a new tumor cell (μm³) */
constexpr real_t kDefaultVolumeNewTumorCell = 2494.0;

/** @brief Default volume of the nucleus of a new tumor cell (μm³) */
constexpr real_t kDefaultVolumeNucleusTumorCell = 540.0;

/** @brief Default fraction of fluid volume in a new tumor cell */
constexpr real_t kDefaultFractionFluidTumorCell = 0.75;

/** @} */ // end of Volume Parameters group


/** @name Volume Relaxation Rates
 *  @brief Volume relaxation rates (min⁻¹) for different cell states
 *  @{
 */

/** @name Alive Cell Volume Relaxation Rates
 *  @brief Relaxation rates for living cells
 *  @{
 */

/** @brief Volume relaxation rate for cytoplasm in alive cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateAliveCytoplasm =0.13/60.;

/** @brief Volume relaxation rate for nucleus in alive cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateAliveNucleus = 0.22/60.;

/** @brief Volume relaxation rate for fluid in alive cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateAliveFluid = 1.3/60.;

/** @} */ // end of Alive Cell Volume Relaxation Rates group

/** @name Necrotic Swelling Volume Relaxation Rates
 *  @brief Relaxation rates for necrotic cells during swelling phase
 *  @{
 */

/** @brief Volume relaxation rate for cytoplasm in necrotic swelling cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticSwelling = 0.0032/60.0;

/** @brief Volume relaxation rate for nucleus in necrotic swelling cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateNucleusNecroticSwelling = 0.013/60.;

/** @brief Volume relaxation rate for fluid in necrotic swelling cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateFluidNecroticSwelling = 0.050/60.0;

/** @} */ // end of Necrotic Swelling Volume Relaxation Rates group

/** @name Necrotic Lysed Volume Relaxation Rates
 *  @brief Relaxation rates for necrotic cells during lysis phase
 *  @{
 */

/** @brief Volume relaxation rate for cytoplasm in necrotic lysed cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticLysed = 0.0032/60.00;

/** @brief Volume relaxation rate for nucleus in necrotic lysed cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateNucleusNecroticLysed = 0.013/60.;

/** @brief Volume relaxation rate for fluid in necrotic lysed cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateFluidNecroticLysed = 0.050/60.0;

/** @} */ // end of Necrotic Lysed Volume Relaxation Rates group
/** @} */ // end of Volume Relaxation Rates group
/** @} */ // end of Tumor Cell Hyperparameters group


/** @name General Simulation Hyperparameters
 *  @brief Core simulation parameters and timing
 *  @{
 */

/** @brief Seed for random number generation */
constexpr int kSeed =3;

/** @name Time Steps
 *  @brief Different time steps for various simulation processes
 *  @{
 */

/** @brief Time step for substance secretion/consumption (minutes) */
constexpr real_t kDtSubstances = 0.01;

/** @brief Time step for cell mechanics (minutes) */
constexpr real_t kDtMechanics = 0.1;

/** @brief Time step for cell cycle processes (minutes) */
constexpr real_t kDtCycle = 6.0;

/** @brief General time step for the simulation (same as kDtMechanics) 
 *  @warning Do not modify this line
 */
constexpr real_t kDt = kDtMechanics;

/** @brief Number of steps per cycle step 
 *  @warning Do not modify this line. Computed to avoid errors with fmod
 */
constexpr int kStepsPerCycle = kDtCycle / kDt;

/** @} */ // end of Time Steps group

/** @name Simulation Duration and Output
 *  @brief Parameters controlling simulation length and output frequency
 *  @{
 */

/** @brief Output summary interval (every 12 hours in simulation time) */
constexpr int kOutputCsvInterval = 12*60/kDt;

/** @brief Total simulation time in minutes (30 days) */
constexpr int kTotalMinutesToSimulate = 30*24*60;

/** @} */ // end of Simulation Duration and Output group

/** @name Spatial Parameters
 *  @brief Parameters defining the simulation space
 *  @{
 */

/** @brief Length of the bounded simulation space (micrometers) */
constexpr int kBoundedSpaceLength = 1000;

/** @brief Initial radius of the spherical tumor (micrometers) */
constexpr real_t kInitialRadiusTumor = 150;

/** @} */ // end of Spatial Parameters group


/** @name Apoptosis Parameters
 *  @brief Parameters for apoptotic cell behavior
 *  @{
 */

/** @brief Volume relaxation rate for cytoplasm in apoptotic cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateCytoplasmApoptotic = 1.0/60.0;

/** @brief Volume relaxation rate for nucleus in apoptotic cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateNucleusApoptotic = 0.35/60.0;

/** @brief Volume relaxation rate for fluid in apoptotic cells (min⁻¹) */
constexpr real_t kVolumeRelaxarionRateFluidApoptotic = 0.0;

/** @brief Time until an apoptotic cell is removed from simulation (minutes) */
constexpr real_t kTimeApoptosis = 8.6*60;

/** @brief Reduction factor for consumption rate of dead cells entering necrosis */
constexpr real_t kReductionConsumptionDeadCells= 0.1;

/** @} */ // end of Apoptosis Parameters group



/** @name Chemical Diffusion Parameters
 *  @brief Parameters for substance diffusion and chemical environment
 *  @{
 */

/** @brief Number of voxels per axis in diffusion grid */
constexpr int kResolutionGridSubstances = 50;

/** @brief Volume of each voxel in the diffusion grid
 *  @warning Do not modify this line
 */
constexpr real_t kVoxelVolume = (kBoundedSpaceLength / kResolutionGridSubstances)*(kBoundedSpaceLength / kResolutionGridSubstances)*(kBoundedSpaceLength/ kResolutionGridSubstances);

/** @name Oxygen Parameters
 *  @brief Parameters specific to oxygen diffusion and behavior
 *  @{
 */

/** @brief Oxygen diffusion coefficient (micrometers²/minute) */
constexpr real_t kDiffusionCoefficientOxygen = 100000;

/** @brief Oxygen decay constant (minutes⁻¹) */
constexpr real_t kDecayConstantOxygen = 0.1;

/** @brief Time step for oxygen diffusion calculations (minutes) */
constexpr real_t kTimeStepOxygen = 0.0005;

/** @brief Reference oxygen level at simulation boundaries (mmHg) */
constexpr real_t kOxygenReferenceLevel = 38.;

/** @brief Initial oxygen concentration in all voxels (mmHg) */
constexpr real_t kInitialOxygenLevel = 38.0;

/** @brief Oxygen saturation level in microenvironment (mmHg) */
constexpr real_t kOxygenSaturation = 30.0;

/** @} */ // end of Oxygen Parameters group

/** @name Immunostimulatory Factor Parameters
 *  @brief Parameters specific to immunostimulatory factor diffusion
 *  @{
 */

/** @brief Immunostimulatory factor diffusion coefficient (micrometers²/minute) */
constexpr real_t kDiffusionCoefficientImmunostimulatoryFactor = 1000;

/** @brief Immunostimulatory factor decay constant (minutes⁻¹) */
constexpr real_t kDecayConstantImmunostimulatoryFactor = 0.016;

/** @brief Time step for immunostimulatory factor diffusion calculations (minutes) */
constexpr real_t kTimeStepImmunostimulatoryFactor = 0.01;

/** @} */ // end of Immunostimulatory Factor Parameters group
/** @} */ // end of Chemical Diffusion Parameters group
/** @name Mechanical Forces Parameters
 *  @brief Parameters controlling cell-cell interaction forces
 *  @{
 */

/** @name Repulsion Forces
 *  @brief Repulsion coefficients between different cell types
 *  @{
 */

/** @brief Repulsion coefficient between tumor cells */
constexpr real_t kRepulsionTumorTumor = 10.0;

/** @brief Repulsion coefficient between CAR-T cells */
constexpr real_t kRepulsionCartCart = 50.0;

/** @brief Repulsion coefficient from CAR-T cells to tumor cells */
constexpr real_t kRepulsionCartTumor = 50.0;

/** @brief Repulsion coefficient from tumor cells to CAR-T cells */
constexpr real_t kRepulsionTumorCart = 10.0;

/** @} */ // end of Repulsion Forces group

/** @name Adhesion Forces
 *  @brief Adhesion coefficients and distance parameters
 *  @{
 */

/** @brief Maximum relative adhesion distance for cell interactions */
constexpr real_t kMaxRelativeAdhesionDistance =1.25;

/** @brief Adhesion coefficient between tumor cells */
constexpr real_t kAdhesionTumorTumor = 0.4;

/** @brief Adhesion coefficient between CAR-T cells */
constexpr real_t kAdhesionCartCart = 0.0;

/** @brief Adhesion coefficient from CAR-T cells to tumor cells */
constexpr real_t kAdhesionCartTumor = 0.0;

/** @brief Adhesion coefficient from tumor cells to CAR-T cells */
constexpr real_t kAdhesionTumorCart = 0.0;

/** @} */ // end of Adhesion Forces group
/** @} */ // end of Mechanical Forces Parameters group

/** @name Computational Parameters
 *  @brief Internal computational parameters and constants
 *  @{
 */

/** @name Adams-Bashforth Coefficients
 *  @brief Coefficients for two-step Adams-Bashforth time derivative approximation
 * 
 *  Position update formula: position(t + dt) ≈ position(t) + dt * [1.5 * velocity(t) - 0.5 * velocity(t - dt)]
 *  @warning Do not change these values
 *  @{
 */

/** @brief Coefficient for current velocity term (dt × 1.5) */
constexpr real_t kDnew= 1.5*kDtMechanics;

/** @brief Coefficient for previous velocity term (dt × -0.5) */
constexpr real_t kDold = -0.5*kDtMechanics;

/** @} */ // end of Adams-Bashforth Coefficients group

/** @brief Length of the mechanics box (micrometers)
 *  @warning Do not change this line
 */
const real_t kLengthBoxMechanics =22;

/** @brief Maximum squared distance for considering cells as neighbors in force calculations (μm²)
 * 
 * Calculated as: (0.1 + cell_radius × kMaxRelativeAdhesionDistance)²
 * Includes 0.1 μm buffer to avoid numerical errors
 * @warning Do not change this line
 */
const real_t kSquaredMaxDistanceNeighborsForce = std::pow(0.1+ std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi) * kMaxRelativeAdhesionDistance,2);

/** @} */ // end of Computational Parameters group
/** @} */ // end of General Simulation Hyperparameters group


/** @name CAR-T Cell Hyperparameters
 *  @brief Parameters controlling CAR-T cell behavior and properties
 *  @{
 */

/** @brief Average maximum time until apoptosis for CAR-T cells (minutes) */
constexpr real_t kAverageMaximumTimeUntillApoptosisCart= kDtCycle* 10.0 * 24.0 * 60.0;

/** @name CAR-T Cell Volume Parameters
 *  @brief Default volume parameters for CAR-T cells
 *  @{
 */

/** @brief Default total volume of a new CAR-T cell (μm³) */
constexpr real_t kDefaultVolumeNewCartCell = 2494.0;

/** @brief Default volume of the nucleus of a new CAR-T cell (μm³) */
constexpr real_t kDefaultVolumeNucleusCartCell = 540.0;

/** @brief Default fraction of fluid volume in a new CAR-T cell */
constexpr real_t kDefaultFractionFluidCartCell = 0.75;

/** @} */ // end of CAR-T Cell Volume Parameters group
/** @} */ // end of CAR-T Cell Hyperparameters group


}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_