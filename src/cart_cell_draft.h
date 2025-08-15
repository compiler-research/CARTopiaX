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

// #ifndef CART_CELL_H_
// #define CART_CELL_H_

// #include "biodynamo.h"
// #include "core/util/log.h"
// #include "core/util/root.h"
// #include "utils_aux.h"
// #include "hyperparams.h"
// #include "tumor_cell.h"

// namespace bdm {


// enum class FasesKi67 : int {// Enum representing the different proliferation phases of the protein Ki67 expression in a cell
//   kKi67Negative = 0,  // Ki67 negative phase
//   kKi67PositivePremitotic = 1,  // Ki67 positive phase before mitosis
//   kKi67PositivePostmitotic = 2,  // Ki67 positive phase just after mitosis, returning to kNegative
// };



// // ─────────────────────────────
// // CartCellState Enum Definition
// // ─────────────────────────────
// /// Enum representing the different states of a tumor cell
// enum class CartCellState : int {
//   kActivated = 0,
//   kApoptotic = 1,
//   kDead = 2
// };


// // ─────────────────────────────
// // CartCell Class Definition
// // ─────────────────────────────
// class CartCell : public Cell {//CHANGE: consumes 1 oxygen per minute //CHANGE: Induces apoptosis include and copy code form tumor_cell entry apoptosis
//   BDM_AGENT_HEADER(CartCell, Cell, 1);
//   public:
//   CartCell() {}
//   explicit CartCell(const real_t3& position) {
//     // Initialize the cell with a given position
//     SetPosition(position);
//     state_ = CartCellState::kActivated;  // Default state for new cells
//     remaining_life_time_ = SamplePositiveGaussian(kAverageLifeTimeCartCell, kStandardDeviationLifeTimeCartCell); //random time following a gaussian distribution
//     timer_last_division_=0;
//     timer_kill_trial_=kTimeKillTrialCart; //new cells are ready to kill
//     exhaustion_level_=0;
//     suppression_level_=0;


//     recognized_antigens_= SampleAntigenPattern(kRecognizedAntigensCart); // Sample a random antigen pattern from the predefined antigen patterns for CAR-T cells
//     this->oxygen_dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid("oxygen"); // Pointer to oxygen diffusion grid
//     this->SetDiameter(kDiameterCartCell);  // Set default diameter
//     this->SetMass(this->GetVolume() * kCartCellDensity);  // Set initial mass
//   }
//   virtual ~CartCell() {}

//   /// Called when a new agent is created (e.g., after cell division)
//   void Initialize(const NewAgentEvent& event) override {
//     Base::Initialize(event);

//     if (auto* mother = dynamic_cast<CartCell*>(event.existing_agent)) {
//       if (event.GetUid() == CellDivisionEvent::kUid) {
//         // Initialize a new CartCell as a daughter of the mother CartCell
//         state_ = CartCellState::kActivated;  // Default state for new cells
//         remaining_life_time_ = SamplePositiveGaussian(kAverageLifeTimeCartCell, kStandardDeviationLifeTimeCartCell); //random time following a gaussian distribution
//         timer_last_division_=0;
//         timer_kill_trial_=mother->timer_kill_trial_; // Inherit the timer for kill trials from the mother cell
//         exhaustion_level_=0;
//         suppression_level_=0;
//         recognized_antigens_= mother->recognized_antigens_;//recognizes the same antigens as the mother cell
//         this->oxygen_dgrid_ = mother->oxygen_dgrid_; // Pointer to oxygen diffusion grid
//         // this->nutrients_dgrid_ = mother->nutrients_dgrid_; // Pointer to nutrients diffusion grid
//         this->SetDiameter(kDiameterCartCell);  // Set default diameter for daughter cell (avoid numerical instabilities)
//         this->SetMass(this->GetVolume() * kCartCellDensity);  // Set initial mass
//         mother->SetMass(mother->GetVolume()*kCartCellDensity); 
//       }
//     }

//   }

//   void SetState(CartCellState state) { state_ = state; }
//   CartCellState GetState() const { return state_; }

//   void SetExpresedAntigens(std::map<std::string, bool> recognized_antigens) { recognized_antigens_ = recognized_antigens; }
//   std::map<std::string, bool> GetExpresedAntigens() const { return recognized_antigens_; }

//   void SetRecognizedAntigens(const std::map<std::string, bool>& recognized_antigens) { recognized_antigens_ = recognized_antigens; }
//   std::map<std::string, bool> GetRecognizedAntigens() const { return recognized_antigens_; }

//   void SetOxygenDGrid(DiffusionGrid* oxygen_dgrid) { oxygen_dgrid_ = oxygen_dgrid; }
//   DiffusionGrid* GetOxygenDGrid() const { return oxygen_dgrid_; }


//   void SetRemainingLifeTime(int remaining_life_time) {remaining_life_time_ = remaining_life_time;}
//   int GetRemainingLifeTime() const { return remaining_life_time_; }

//   void SetTimerLastDivision(int timer_last_division) { timer_last_division_ = timer_last_division; }
//   int GetTimerLastDivision() const { return timer_last_division_; }

//   void SetTimerKillTrial(int timer_kill_trial) { timer_kill_trial_ = timer_kill_trial; }
//   int GetTimerKillTrial() const { return timer_kill_trial_; }

//   void SetTimerState(int timer_state) { timer_state_ = timer_state; }
//   int GetTimerState() const { return timer_state_; }

//   void SetExhaustionLevel(real_t exhaustion_level) { exhaustion_level_ = exhaustion_level; }
//   real_t GetExhaustionLevel() const { return exhaustion_level_; }

//   void SetSuppressionLevel(real_t suppression_level) { suppression_level_ = suppression_level; }
//   real_t GetSuppressionLevel() const { return suppression_level_; }

//   void SetTimerUnderSuppression(int timer_under_suppression) { timer_under_suppression_ = timer_under_suppression; }
//   int GetTimerUnderSuppression() const { return timer_under_suppression_; }


//   /// Returns the diffusion grid for oxygen
//   DiffusionGrid* GetOxygenDiffusionGrid() const { return oxygen_dgrid_; }
//   /// Returns the diffusion grid for nutrients
//   // DiffusionGrid* GetNutrientsDiffusionGrid() const { return nutrients_dgrid_; }


//  private:
//   CartCellState state_;
//   int remaining_life_time_;//remaining steps life time of the cell. FIrst used to turn the cell apoptotic, then to remove it from the simulation
//   int timer_last_division_;  // Timer since the last division
//   int timer_kill_trial_;  // Timer sice the last kill trial
//   int timer_state_; // Timer for the current state (used for apoptosis and dead states)
//   real_t exhaustion_level_;  // CAR-T exhaustion level (Permanent damage)
//   real_t suppression_level_;  // CAR-T suppression level (Cyotokines temporal damage)
//   int timer_under_suppression_;  // Counter of steps under high suppression
//   std::map<std::string, bool> recognized_antigens_; // Recognized antigens
//   DiffusionGrid* oxygen_dgrid_;  // Pointer to the oxygen diffusion grid
//   // DiffusionGrid* nutrients_dgrid_;  // Pointer to the nutrients diffusion grid
//   //CHANGE add cytokines grids
// };

// // ─────────────────────────────
// // Behavior: StateControlGrowProliferate
// // ─────────────────────────────
// struct StateControlAndCartDivision : public Behavior {
//   BDM_BEHAVIOR_HEADER(StateControlAndCartDivision, Behavior, 1);

//   StateControlAndCartDivision() { AlwaysCopyToNew(); }
//   virtual ~StateControlAndCartDivision() {}

//   /// Main behavior executed at each simulation step
//   void Run(Agent* agent) override {
//     auto* sim = Simulation::GetActive();
//     auto* random = sim->GetRandom();

//     if (auto* cell = dynamic_cast<CartCell*>(agent)) {

//       // Nutrients and oxygen levels
//       real_t3 current_position = cell->GetPosition();
//       auto* oxygen_dgrid = cell->GetOxygenDiffusionGrid();  // Pointer to the oxygen diffusion grid
//       // auto* nutrients_dgrid = cell->GetNutrientsDiffusionGrid();  // Pointer to the nutrients diffusion grid
//       real_t oxygen_level = oxygen_dgrid->GetValue(current_position);
//       // real_t nutrients_level = nutrients_dgrid->GetValue(current_position);
//       // Get the cytokines levels CHANGE

      
//       //compute Suppression level
//       real_t suppression_level =0.;// f(cytokines) //CHANGE
//       //Increase suppression level based on oxygen and nutrients levels
//       if (oxygen_level < kThresholdOxygenLevelSuppression) { suppression_level += 0.1;}
//       // if (nutrients_level < kThresholdNutrientsLevelSuppression) { suppression_level += 0.1; }
//       cell->SetSuppressionLevel(suppression_level);

//       //Exaustion level from suppression level
//       if (suppression_level > kThresholdHighSuppression) {
//         cell->SetTimerUnderSuppression(cell->GetTimerUnderSuppression() + 1); // Increment timer under suppression
//         if (cell->GetTimerUnderSuppression() >= kThresholdStepsUnderHighSuppression) { // If under suppression for too long
//           real_t previous_exhaustion_level = cell->GetExhaustionLevel();
//           real_t diff = suppression_level- kThresholdHighSuppression; // Difference from the threshold
//           cell->SetExhaustionLevel(std::min(1.0, previous_exhaustion_level+ SamplePositiveGaussian(diff*diff,0.005))); // Increase exhaustion level based on suppression level
//         }
//       } else {
//         cell->SetTimerUnderSuppression(0); // Reset timer under suppression
//       }
//       //CHANGE increase exhaustion level based on certain cytokines that increase it directly

//       switch (cell->GetState())
//       {
//       case CartCellState::kActivated: // Cell is in nomral activated state
//         cell->SetTimerLastDivision(cell->GetTimerLastDivision() + 1);  // Increment timer_last_division
//         if (cell->GetRemainingLifeTime() <= 0) {
//           cell->SetState(CartCellState::kApoptotic);  // Transition to apoptotic state
//           cell->SetTimerState(0);  // Reset timer_state
//         } else {
//           // Decrease remaining life time
//           auto aux_exp=std::exp(3*cell->GetExhaustionLevel());
//           real_t decrement= SamplePositiveGaussian(aux_exp, 0.2*aux_exp); // Decrease life time based on exhaustion level
//           cell->SetRemainingLifeTime(cell->GetRemainingLifeTime() - decrement);
//           // Check if the cell can divide
//           if (cell->GetTimerLastDivision() >= kTimeLastDivisionCartCell) {
//             if (random->Uniform(0.0, 1.0) < kBaseProbabilityDivideCartCell) {//CHANGE: probability needs to depend on the nutrients and cytokines levels
//               cell->SetTimerLastDivision(0);  // Reset timer_last_division
//               // Create a new daughter cell
//               cell->SetDiameter(std::cbrt(2.0) *kDiameterCartCell);  // adjust the diameter so that the volume of the two new cells is conserved is the same as before division
//               cell->Divide();  // Perform cell division
//               cell->SetDiameter(kDiameterCartCell);  // Reset diameter to original value to avoid insteabilities
//             }
//           }
//         }
//         break;
//       case CartCellState::kApoptotic:
//           // Apoptosis induced cells die after a certain time
//           cell->SetTimerState(cell->GetTimerState() + 1);  // Increase timer_state to track duration
//           if (cell->GetTimerState() >= kTimeApoptosisInducedCart ) {
//             cell->SetState(CartCellState::kDead);
//             cell->SetTimerState(0);  // Reset timer_state
//           }
//           break;
//       case CartCellState::kDead:
//           cell->SetTimerState(cell->GetTimerState() + 1);  // Increase timer_state to track duration
//           if (cell->GetTimerState() >= kTimeDeadCart ) {
//             //remove the cell from the simulation
//             auto* ctxt = sim->GetExecutionContext();
//             ctxt->RemoveAgent(agent->GetUid());
//           }
//           break;
//       default:
//         Log::Error("StateControlGrowProliferate::Run", "Unknown CartCellState");
//         break;
//       }
//     } else {
//       Log::Error("StateControlGrowProliferate::Run", "SimObject is not a CartCell");
//     }
//   }
// };

// struct CartMigration : public Behavior{
//   BDM_BEHAVIOR_HEADER(CartMigration, Behavior, 1);

//   CartMigration() { AlwaysCopyToNew(); }
//   virtual ~CartMigration() {}

//   /// Main behavior executed at each simulation step
//   void Run(Agent* agent) override {
//     auto* sim = Simulation::GetActive();
//     auto* random = sim->GetRandom();
//     if (auto* cell = dynamic_cast<CartCell*>(agent)) {

//       real_t3 current_position = cell->GetPosition();
      
//       // auto* oxygen_dgrid = cell->GetOxygenDiffusionGrid(); //CHANGE: get cytokines grids
//       real_t speed = kMaxSpeedCartCell;  // Default speed

//       // Adjust speed based on suppression and exhaustion levels
//       // CHANGE this equation
//       speed = kMaxSpeedCartCell * (0.5*(1-cell->GetSuppressionLevel()) + 0.5*(1-cell->GetExhaustionLevel()));


//       Real3 direction;
//       switch (cell->GetState())
//       {
//       case CartCellState::kActivated:
//         //the speed should is reduced in hypoxic conditions because of the lack of oxygen in the equation used prevously
//         // dcytokines->GetGradient(current_position, &direction);  // returns normalized gradient towards the cytokines source CHANGE
//         direction = {// Move randomly, CHANGE: this should be a gradient towards the cytokines source
//             random->Uniform(-0.1, 0.1),
//             random->Uniform(-0.1, 0.1),
//             random->Uniform(-0.1, 0.1)};;
//         break;
//       case CartCellState::kApoptotic:
//         // Apoptosis induced cells may not move or move less
//         speed = speed * 0.2;  // Further reduced speed for apoptosis induced cells
//         direction = {
//             random->Uniform(-0.1, 0.1),
//             random->Uniform(-0.1, 0.1),
//             random->Uniform(-0.1, 0.1)};
//         break;
//       case CartCellState::kDead:
//         // Dead cells do not move
//         return;
//       default:
//         Log::Error("CartMigration::Run", "Unknown CartCellState");
//         return;
//       }

//       // Update the cell's position
//       cell->SetTractorForce(direction * speed);
//     }
//   }
// };

// // struct CartReleaseChemicals : public Behavior {
// //   BDM_BEHAVIOR_HEADER(CartReleaseChemicals, Behavior, 1);

// //   CartReleaseChemicals() { AlwaysCopyToNew(); }
// //   virtual ~CartReleaseChemicals() {}

// //   /// Main behavior executed at each simulation step
// //   void Run(Agent* agent) override {
// //     auto* sim = Simulation::GetActive();
// //     auto* random = sim->GetRandom();
// //     //nutrients= sim->GetEnvironment()->GetNutrients(); //CHANGE: Get the current nutrients level from the environment and effect the chemicals released
// //     if (auto* cell = dynamic_cast<CartCell*>(agent)) {

// //       switch (cell->GetState())
// //       {
// //       case CartCellState::kProliferative:

// //         break;
// //       case CartCellState::kGrowing:

// //         break;
// //       case CartCellState::kHypoxic:

// //         break;
// //       case CartCellState::kApoptosisInduced:

// //         break;
// //       case CartCellState::kDead:
// //         return;
// //       default:
// //         Log::Error("CartReleaseChemicals::Run", "Unknown CartCellState");
// //         return;
// //       }
// //     }
// //   }
// // };

// class CartAttackBehavior : public Behavior {
//   BDM_BEHAVIOR_HEADER(CartAttackBehavior, Behavior, 1);

//  public:
//   CartAttackBehavior() { AlwaysCopyToNew(); }
//   virtual ~CartAttackBehavior() {}

//   void Run(Agent* agent) override {
//     if (auto* cart = dynamic_cast<CartCell*>(agent)) {
//       cart->SetTimerKillTrial(cart->GetTimerKillTrial() + 1); // Increase the timer for kill trials
//       if (cart->GetState() == CartCellState::kActivated && cart->GetTimerKillTrial() >= kTimeKillTrialCart){// If the cell is activated and ready to attack
//         auto* sim= Simulation::GetActive();
//         auto* random = sim->GetRandom();

//         auto* ctxt = sim->GetExecutionContext();
//         const auto& cart_pos = cart->GetPosition();
//         real_t radius = (cart->GetDiameter() + kMaxDiameterTumorCell) * 0.5;//contact distance with tumor cells
//         real_t squared_radius = radius * radius + 0.01; // Add a small value to avoid numerical issues

//         TumorCell* closest_target = nullptr;
//         real_t min_sq_dist = std::numeric_limits<real_t>::max();

//         // Lambda to find the closest alive tumor cell within contact distance
//         auto lambda_neighbors = [&](Agent* neighbor, real_t squared_distance) {
//           if (auto* tumor = dynamic_cast<TumorCell*>(neighbor)) {
//             real_t contact_dist = (cart->GetDiameter() + tumor->GetDiameter()) * 0.5;
//             // Find the closest target within contact distance(+ epsilon) and exclude dead tumor cells (apoptotic are included as possible targets)
//             auto recognized_antigens = cart->GetRecognizedAntigens();
//             // auto tumor_antigens = tumor->GetExpressedAntigens();

//             bool recognized_as_target = false;

//             // Check if CART recognizes the tumor cell's antigens to target it
//             for (const auto& pair : recognized_antigens) {
//               const std::string& antigen = pair.first;
//               bool is_recognized = pair.second;
//               if (is_recognized) {
//                 auto it = tumor_antigens.find(antigen);
//                 real_t prob = (it != tumor_antigens.end()) ? it->second : 0.0;
//                 if (random->Uniform(0.0, 1.0) < prob) {
//                   recognized_as_target = true;//if it detects at any of the antigens, the tumor cell is recognized as a target
//                   break;
//                 }
//               }
//             }

//             if (recognized_as_target && // If the tumor cell is recognized as a target to be killed
//                 squared_distance <= contact_dist * contact_dist + 0.01 && // Ensure both cells are within contact distance( + epsilon for numerical stability)
//                 squared_distance < min_sq_dist //&&// Looking for the closest target
//                 // tumor->GetState() != TumorCellState::kDead
//               ) {// the target must not be dead (it can be apoptotic though)
//               closest_target = tumor;
//               min_sq_dist = squared_distance;
//             }
//           }
//         };

//         // sim->GetEnvironment()->ForcedUpdate();

//         // auto functor = bdm::L2F(lambda_neighbors);
//         // ctxt->ForEachNeighbor(functor, cart_pos, squared_radius);

//         // if (closest_target) {//there is a target to attack
//         //   cart->SetTimerKillTrial(0); // Reset the timer for kill trials after an attack
//         //   // Attack the closest target
          
//         //   real_t prob = (1-closest_target->GetBaseImmunogenicity())*(0.7*(1-cart->GetSuppressionLevel())+0.3*(1-cart->GetExhaustionLevel())); //CHANGE: make a kill probability dependent of suppression and exhaustion levels and tumor immunogenicity
//         //   if (random->Uniform(0.0, 1.0) < prob) {
//         //     closest_target->SetState(TumorCellState::kApoptosisInduced);
//         //     closest_target->SetTimerState(0); // Reset the timer for the target cell
//         //   }
//         // }
//       }
//     }
//   }
// };


// }  // namespace bdm

// #endif  // CART_CELL_H_
