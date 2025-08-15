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

#ifndef TUMOR_HYPERPARAMS_H_
#define TUMOR_HYPERPARAMS_H_

#include <map>
#include <string>

namespace bdm {



// ─────────────────────────────
// TumorCell Hyperparameters
// ─────────────────────────────

constexpr real_t kRateSecretionImmunostimulatoryFactor= 10.0; // Rate of secretion of immunostimulatory factor of tumor cells per minute
constexpr real_t kSaturationDensityImmunostimulatoryFactor = 1.0; // Saturation density of immunostimulatory factor for tumor cells
constexpr real_t kOncoproteinMean = 1.0; // Mean level of oncoprotein expression in tumor cells
constexpr real_t kOncoproteinStandardDeviation = 0.25; // Standard deviation of oncoprotein expression in tumor cells
constexpr real_t kOxygenSaturationInProliferation = 38.0; // Oxygen saturation level in tumor cells for proliferation
constexpr real_t kOxygenLimitForProliferation = 10.0; // Limit of oxygen level for tumor cell proliferation
constexpr real_t kOxygenLimitForNecrosis = 5.0; // Limit of oxygen to start causing necrosis
constexpr real_t kOxygenLimitForNecrosisMaximum= 2.5; // Limit of oxygen to maximum necrosis probability
// constexpr real_t kTransitionRateKi67[] = {// Transition rates for Ki67 expression states in min**-1
//     1.0 / (3.62 * 60.0),  // Rate from 0 to 1 (ki67 negative to ki67 positive pre-mitotic)
//     1.0 / (13.0 * 60.0),  // Rate from 1 to 2 (ki67 positive pre-mitotic to ki67 positive post-mitotic)
//     1.0 / (2.5 * 60.0)    // Rate from 2 to 0 (ki67 positive post-mitotic to ki67 negative)
// };
constexpr real_t kTimeLysis = 60*24*60.; // Time in minutes until a lysed necrotic cell is removed from the simulation
constexpr real_t kDivisionRate = 0.02717 / 60.0; // Rate of cell division in min**-1
constexpr real_t kMaximumNecrosisRate= 1.0 / (6.0 * 60.0); // Maximum rate per minute of necrosis for tumor cells in case of hypoxia with 0 oxygen
constexpr real_t kDefaultOxygenConsumption = 10.0; // Default oxygen consumption rate of tumor cell
//Volume parameters
constexpr real_t kDefaultVolumeNewTumorCell = 2494.0; // Default total volume of a new tumor cell in μm³
constexpr real_t kDefaultVolumeNucleusTumorCell = 540.0; // Default volume of the nucleus of a new tumor cell in μm³
constexpr real_t kDefaultFractionFluidTumorCell = 0.75; // Default fraction of fluid volume in a new tumor cell


//volume relaxation rate (min^-1) for each state
constexpr real_t kVolumeRelaxarionRateAliveCytoplasm =0.13/60.;// 0.27/ 60.0;
constexpr real_t kVolumeRelaxarionRateAliveNucleus = 0.22/60.;//0.33/60.
constexpr real_t kVolumeRelaxarionRateAliveFluid = 1.3/60.;//3.0/60.

constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticSwelling = 0.0032/60.0;
constexpr real_t kVolumeRelaxarionRateNucleusNecroticSwelling = 0.013/60.;
constexpr real_t kVolumeRelaxarionRateFluidNecroticSwelling = 0.050/60.0;

constexpr real_t kVolumeRelaxarionRateCytoplasmNecroticLysed = 0.0032/60.00;
constexpr real_t kVolumeRelaxarionRateNucleusNecroticLysed = 0.013/60.;
constexpr real_t kVolumeRelaxarionRateFluidNecroticLysed = 0.050/60.0;


// ─────────────────────────────
// General Hyperparameters
// ─────────────────────────────

constexpr int kSeed =3; // Seed for random number generation

constexpr real_t kDtSubstances = 0.01; // 0.01 minutes time step for substances secretion/consumption
constexpr real_t kDtMechanics = 0.1; // 0.1 minutes time step for the cell mechanics
constexpr real_t kDtCycle = 6.0; // 6 minutes time step for the cell cycle

constexpr real_t kDt = kDtMechanics; // General time step for the simulation: it is the same as kDtMechanics, do not modify this line
constexpr int kStepsPerCycle = kDtCycle / kDt; // Number of steps per cycle step, do not modify this line. Needs to be computed to avoid errors with fmod

constexpr int kOutputCsvInterval = 12*60/kDt;// Output little summary each half a day


constexpr int kTotalMinutesToSimulate = 30*24*60; //30 * 24 * 60; // Total simulation time in minutes (30 days)
constexpr int kBoundedSpaceLength = 1000; // Length of the bounded space in micrometers
constexpr real_t kInitialRadiusTumor = 150; // Initial radius of the spherical tumor (group of cancer cells) in micrometers


constexpr real_t kVolumeRelaxarionRateCytoplasmApoptotic = 1.0/60.0;
constexpr real_t kVolumeRelaxarionRateNucleusApoptotic = 0.35/60.0;
constexpr real_t kVolumeRelaxarionRateFluidApoptotic = 0.0;
constexpr real_t kTimeApoptosis = 8.6*60; // Time in minutes until an apoptotic cell is removed from the simulation
constexpr real_t kReductionConsumptionDeadCells= 0.1; // Reduction of consumption rate of dead cells when they enter necrosis



//Chemicals
constexpr int kResolutionGridSubstances = 50; //50 // voxels per axis
constexpr real_t kVoxelVolume = (kBoundedSpaceLength / kResolutionGridSubstances)*(kBoundedSpaceLength / kResolutionGridSubstances)*(kBoundedSpaceLength/ kResolutionGridSubstances); //Do not modify this line
constexpr real_t kDiffusionCoefficientOxygen = 100000; // 100000 micrometers^2/minute
constexpr real_t kDecayConstantOxygen = 0.1; // 0.1 minutes^-1
constexpr real_t kDiffusionCoefficientImmunostimulatoryFactor = 1000; // 1000 micrometers^2/minute
constexpr real_t kDecayConstantImmunostimulatoryFactor = 0.016; // 0.016 minutes^-1
constexpr real_t kTimeStepOxygen = 0.0005; // 0.001 minutes CHANGE
constexpr real_t kTimeStepImmunostimulatoryFactor = 0.01; // 0.01 minutes
constexpr real_t kOxygenReferenceLevel = 38.; // Reference level of oxygen at the boundaries of the simulation space in mmHg
constexpr real_t kInitialOxygenLevel = 38.0; // Initial voxel concentration of oxygen in mmHg
constexpr real_t kOxygenSaturation = 30.0; //30.0 // Oxygen saturation in mmHg in microenvironment
//Forces
constexpr real_t kRepulsionTumorTumor = 10.0; // Repulsion coeficient between tumor cells
constexpr real_t kRepulsionCartCart = 50.0;   // Repulsion coeficient between CAR-T cells
constexpr real_t kRepulsionCartTumor = 50.0;  // Repulsion coeficient between CAR-T cells and tumor cells
constexpr real_t kRepulsionTumorCart = 10.0;  // Repulsion coeficient between tumor cells and CAR-T cells
constexpr real_t kMaxRelativeAdhesionDistance =1.25; // Maximum relative adhesion distance for cell interactions
constexpr real_t kAdhesionTumorTumor = 0.4; // Adhesion coeficient between tumor cells
constexpr real_t kAdhesionCartCart = 0.0;   // Adhesion coeficient between CAR-T cells
constexpr real_t kAdhesionCartTumor = 0.0;  // Adhesion coeficient between CAR-T cells and tumor cells
constexpr real_t kAdhesionTumorCart = 0.0;  // Adhesion coeficient between tumor cells and CAR-T cells

//Do not change
//coefficientes for the two step Adams-Bashforth approximation of the time derivative for position 
//position(t + dt) ≈ position(t) + dt * [ 1.5 * velocity(t) - 0.5 * velocity(t - dt) ]
constexpr real_t kDnew= 1.5*kDtMechanics; //dt*1.5
constexpr real_t kDold = -0.5*kDtMechanics; // dt*(-0.5)

//Do not change this line
const real_t kLengthBoxMechanics =22; // Length of the box for mechanics in micrometers

//Max Distance for considering two cells as neighbours for force calculations in μm
//Do not change this line
const real_t kSquaredMaxDistanceNeighborsForce = std::pow(0.1+ std::cbrt(kDefaultVolumeNewTumorCell * 6 / Math::kPi) * kMaxRelativeAdhesionDistance,2);// (twice biggest cell radius (in case to cells tha maximum size encounter each other) times kMaxRelativeAdhesionDistance + 0.1 to avoid mismatch because of numerical errors)**2


// ─────────────────────────────
// CAR-T Cell Hyperparameters
// ───────────────────────────── 
constexpr real_t kAverageMaximumTimeUntillApoptosisCart= kDtCycle* 10.0 * 24.0 * 60.0;
//Volume parameters
constexpr real_t kDefaultVolumeNewCartCell = 2494.0; // Default total volume of a new CAR-T cell in μm³
constexpr real_t kDefaultVolumeNucleusCartCell = 540.0; // Default volume of the nucleus of a new CAR-T cell in μm³
constexpr real_t kDefaultFractionFluidCartCell = 0.75; // Default fraction of fluid volume in a new CAR-T cell


// ─────────────────────────────
// CAR-T Cell Hyperparameters
// ───────────────────────────── 


// //////////////////////////////////////////////////////old
// constexpr real_t kDiameterCartCell = 10.0;  // Diameter in μm
// constexpr real_t kVolumeGrowthCartCell = 800.0;         // Volume increment per step in μm³
// constexpr real_t kBaseProbabilityDivideCartCell = 0.05;      // Probability of cell division per step
// constexpr int kTimeApoptosisInducedCart = 30;       // Time until apoptosis-induced cell dies in steps
// constexpr int kTimeDeadCart = 40;                   // Time until dead cell is removed in steps
// constexpr real_t kTimeLastDivisionCartCell = 250; // Time until the cell can divide again after division in steps
// constexpr real_t kTimeKillTrialCart = 100;          // Time until the cell can perform a kill trial again in steps
// constexpr real_t kMaxSpeedCartCell = 5.0;             // Random displacement motility scale in μm
// constexpr real_t kCartCellDensity = 1.0; // Initial mass of a new tumor cell
// constexpr real_t kAverageLifeTimeCartCell = 500; // Average life time of a CAR-T cell in steps
// constexpr real_t kStandardDeviationLifeTimeCartCell = 20; // Standard deviation of the life time of a CAR-T cell in steps
// constexpr real_t kThresholdHighSuppression = 0.5; // Threshold for CAR-T cell tu be considered under high suppression
// constexpr real_t kThresholdStepsUnderHighSuppression = 10; // Minimum amount of steps a CAR-T cell has to be under high suppression for exhaustion to be increased
// constexpr real_t kThresholdOxygenLevelSuppression = 0.2; // Threshold for hypoxia in CAR-T cells
// constexpr real_t kThresholdNutrientsLevelSuppression = 0.05; // Threshold for low nutrients in CAR-T cells
// inline std::vector<std::pair<std::map<std::string, bool>, float>> kRecognizedAntigensCart = {//
//     { 
//       {// Antigen recognition this type of for CAR-T cells; 
//       {"HER2", true},
//       {"MUC1", false},
//       {"EpCAM", false},
//       }, 
//       0.6f //proportion of cells that with the previous pattern of antigen recognition
//     },
//     {    
//       {
//       {"HER2", false},
//       {"MUC1", false},
//       {"EpCAM", true},
//       }, 
//       0.4f//40% of cells recognize only EpCAM antigen
//     },//the proportions need to add up to 1.0
// };



}  // namespace bdm

#endif  // TUMOR_HYPERPARAMS_H_