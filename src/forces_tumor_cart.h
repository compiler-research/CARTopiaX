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

#ifndef FORCES_TUMOR_CART_H_
#define FORCES_TUMOR_CART_H_
#include "core/interaction_force.h"
#include "core/operation/mechanical_forces_op.h"
#include "biodynamo.h"
#include "core/util/log.h"
#include "core/util/root.h"
#include "hyperparams.h"
#include "utils_aux.h"
#include "tumor_cell.h"

namespace bdm {

class InteractionVelocity : public InteractionForce {
 public:
  InteractionVelocity() = default;
  ~InteractionVelocity() override = default;

  Real4 Calculate(const Agent* lhs, const Agent* rhs) const override;

  InteractionForce* NewCopy() const override;
};

}  // namespace bdm

#endif  // FORCES_TUMOR_CART_H_