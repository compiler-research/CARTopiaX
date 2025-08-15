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

#ifndef TUMOR_HYPERPARAMS_H_
#define TUMOR_HYPERPARAMS_H_

#include "core/real_t.h"
#include "core/util/math.h"
#include <cmath>
#include <cstddef>

namespace bdm {

/// This file contains hyperparameters used in the simulation. Change: In a
/// future version it needs to be changed into a params file with no need to be
/// recompiled

///
/// TumorCell Hyperparameters
///

/// Rate of secretion of immunostimulatory factor of tumor cells per minute
constexpr real_t kRateSecretionImmunostimulatoryFactor = 10.0;
/// Saturation density of immunostimulatory factor for tumor cells
constexpr real_t kSaturationDensityImmunostimulatoryFactor = 1.0;
/// Mean level of oncoprotein expression in tumor cells
constexpr real_t kOncoproteinMean = 1.0;
/// Standard deviation of oncoprotein expression in tumor cells
constexpr real_t kOncoproteinStandardDeviation = 0.25;
/// Oxygen saturation level in tumor cells for proliferation
constexpr real_t kOxygenSaturationInProliferation = 38.0;
/// Limit of oxygen level for tumor cell proliferation
constexpr real_t kOxygenLimitForProliferation = 10.0;
/// Limit of oxygen to start causing necrosis
constexpr real_t kOxygenLimitForNecrosis = 5.0;
/// Limit of oxygen to maximum necrosis probability
constexpr real_t kOxygenLimitForNecrosisMaximum = 2.5;
/// Time in minutes until a lysed necrotic cell is removed from the simulation
constexpr real_t kTimeLysis = 60 * 24 * 60.;
/// Rate of cell division in min**-1
constexpr real_t kDivisionRate = 0.02717 / 60.0;
/// Maximum rate per minute of necrosis for tumor cells in case of hypoxia with
/// 0 oxygen
constexpr real_t kMaximumNecrosisRate = 1.0 / (6.0 * 60.0);
/// Default oxygen consumption rate of tumor cell
constexpr real_t kDefaultOxygenConsumption = 10.0;
/// Volume parameters
///  Default total volume of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNewTumorCell = 2494.0;
/// Default volume of the nucleus of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNucleusTumorCell = 540.0;
/// Default fraction of fluid volume in a new tumor cell
constexpr real_t kDefaultFractionFluidTumorCell = 0.75;

/// Average time for transformation Random Rate in hours
constexpr real_t kAverageTimeTransformationRandomRate = 38.6;

/// Standard Deviation for transformation Random Rate in hours
constexpr real_t kStandardDeviationTransformationRandomRate = 3.7;

/// volume relaxation rate (min^-1) for each state
constexpr real_t kVolumeRelaxarionRateAliveCytoplasm = 0.13 / 60.;
constexpr real_t kVolumeRelaxarionRateAliveNucleus = 0.22 / 60.;
constexpr real_t kVolumeRelaxarionRateAliveFluid = 1.3 / 60.;

constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticSwelling = 0.0032 / 60.0;
constexpr real_t kVolumeRelaxarionRateNucleusNecroticSwelling = 0.013 / 60.;
constexpr real_t kVolumeRelaxarionRateFluidNecroticSwelling = 0.050 / 60.0;

constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticLysed = 0.0032 / 60.00;
constexpr real_t kVolumeRelaxarionRateNucleusNecroticLysed = 0.013 / 60.;
constexpr real_t kVolumeRelaxarionRateFluidNecroticLysed = 0.050 / 60.0;

/// Thresholds in oncoprotein levels for differentiating 4 cancer cell types
constexpr real_t kThresholdCancerCellType1 = 1.5;
constexpr real_t kThresholdCancerCellType2 = 1.0;
constexpr real_t kThresholdCancerCellType3 = 0.5;
constexpr real_t kThresholdCancerCellType4 = 0.0;

///
/// General Hyperparameters
///
/// Seed for random number generation
constexpr int kSeed = 3;

/// 0.01 minutes time step for substances secretion/consumption
constexpr real_t kDtSubstances = 0.01;
/// 0.1 minutes time step for the cell mechanics
constexpr real_t kDtMechanics = 0.1;
/// 6 minutes time step for the cell cycle
constexpr real_t kDtCycle = 6.0;

/// General time step for the simulation: it is the same as kDtMechanics, do not
/// modify this line
constexpr real_t kDt = kDtMechanics;
/// Number of steps per cycle step, do not modify this line. Needs to be
/// computed to avoid errors with fmod
constexpr int kStepsPerCycle = kDtCycle / kDt;

/// Epsilon for avoiding division by 0
constexpr real_t kEpsilon = 1e-10;

/// kEpsilon distance
constexpr real_t kEpsilonDistance = 1e-5;

/// Output little summary each half a day
constexpr int kOutputCsvInterval = 12 * 60 / kDt;

/// Total simulation time in minutes (30 days)
constexpr int kTotalMinutesToSimulate = 30 * 24 * 60;
/// Length of the bounded space in micrometers
constexpr int kBoundedSpaceLength = 1000;
/// Initial radius of the spherical tumor (group of cancer cells) in micrometers
constexpr real_t kInitialRadiusTumor = 150;

constexpr real_t kVolumeRelaxarionRateCytoplasmApoptotic = 1.0 / 60.0;
constexpr real_t kVolumeRelaxarionRateNucleusApoptotic = 0.35 / 60.0;
constexpr real_t kVolumeRelaxarionRateFluidApoptotic = 0.0;
/// Time in minutes until an apoptotic cell is removed from the simulation
constexpr real_t kTimeApoptosis = 8.6 * 60;
/// Reduction of consumption rate of dead cells when they enter necrosis
constexpr real_t kReductionConsumptionDeadCells = 0.1;

/// Constant for division by 2.0
constexpr real_t kHalf = 2.0;

/// Chemicals
///  Number of voxels per axis for the substances grid
constexpr int kResolutionGridSubstances = 50;
/// Volume of a single voxel in μm³ (do not modify this line)
constexpr real_t kVoxelVolume =
    (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances) *
    (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances) *
    (static_cast<real_t>(kBoundedSpaceLength) / kResolutionGridSubstances);
/// Diffusion coefficient of oxygen in μm²/min
constexpr real_t kDiffusionCoefficientOxygen = 100000;
/// Decay constant of oxygen in min⁻¹
constexpr real_t kDecayConstantOxygen = 0.1;
/// Diffusion coefficient of immunostimulatory factor in μm²/min
constexpr real_t kDiffusionCoefficientImmunostimulatoryFactor = 1000;
/// Decay constant of immunostimulatory factor in min⁻¹
constexpr real_t kDecayConstantImmunostimulatoryFactor = 0.016;
/// Reference level of oxygen at the boundaries in mmHg
constexpr real_t kOxygenReferenceLevel = 38.;
/// Initial oxygen concentration in each voxel in mmHg
constexpr real_t kInitialOxygenLevel = 38.0;
/// Oxygen saturation in the microenvironment in mmHg
constexpr real_t kOxygenSaturation = 30.0;
/// Forces
///  Repulsion coeficient between tumor cells
constexpr real_t kRepulsionTumorTumor = 10.0;
/// Repulsion coeficient between CAR-T cells
constexpr real_t kRepulsionCartCart = 50.0;
/// Repulsion coeficient between CAR-T cells and tumor cells
constexpr real_t kRepulsionCartTumor = 50.0;
/// Repulsion coeficient between tumor cells and CAR-T cells
constexpr real_t kRepulsionTumorCart = 10.0;
/// Maximum relative adhesion distance for cell interactions
constexpr real_t kMaxRelativeAdhesionDistance = 1.25;
/// Adhesion coeficient between tumor cells
constexpr real_t kAdhesionTumorTumor = 0.4;
/// Adhesion coeficient between CAR-T cells
constexpr real_t kAdhesionCartCart = 0.0;
/// Adhesion coeficient between CAR-T cells and tumor cells
constexpr real_t kAdhesionCartTumor = 0.0;
/// Adhesion coeficient between tumor cells and CAR-T cells
constexpr real_t kAdhesionTumorCart = 0.0;

/// Do not change
// coefficientes for the two step Adams-Bashforth approximation of the time
// derivative for position position(t + dt) ≈ position(t) + dt * [ 1.5 *
// velocity(t) - 0.5 * velocity(t - dt) ]
/// Coefficient for the current time step in the Adams-Bashforth method (dt
/// * 1.5)
constexpr real_t kDnew = 1.5 * kDtMechanics;
/// Coefficient for the previous time step in the Adams-Bashforth method (dt *
/// -0.5)
constexpr real_t kDold = -0.5 * kDtMechanics;

/// Large time to avoid division by 0
constexpr real_t kTimeTooLarge = 1e100;

/// Minutes in an hour
constexpr real_t kMinutesInAnHour = 60.0;
/// Hours in a day
constexpr real_t kHoursInADay = 24.0;

/// Do not change this line
const size_t gKLengthBoxMechanics = 22;

/// Max Distance for considering two cells as neighbours for force calculations
/// in μm Do not change this line
/// (twice biggest cell radius (in case to cells tha maximum size
/// encounter each other) times kMaxRelativeAdhesionDistance + 0.1 to
/// avoid mismatch because of numerical errors)**2
const real_t gKSquaredMaxDistanceNeighborsForce =
    std::pow(0.1 + std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi) *
                       kMaxRelativeAdhesionDistance,
             2);

///
/// CAR-T Cell Hyperparameters
///
constexpr real_t kAverageMaximumTimeUntillApoptosisCart =
    kDtCycle * 10.0 * 24.0 * 60.0;
/// Volume parameters
///  Default total volume of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNewCarTCell = 2494.0;
/// Default volume of the nucleus of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNucleusCarTCell = 540.0;
/// Default fraction of fluid volume in a new CAR-T cell
constexpr real_t kDefaultFractionFluidCarTCell = 0.75;

}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_
