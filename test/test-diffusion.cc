#include "diffusion/diffusion_thomas_algorithm.h"
#include "params/hyperparams.h"
#include <gtest/gtest.h>
#include "biodynamo.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {
class DiffusionThomasAlgorithmTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Register before constructing the Simulation, as Simulate() does.
    // BioDynaMo materialises the parameter groups while the Simulation is
    // being built, so a group registered afterwards is never visible to
    // Param::Get<SimParam>() -- which then returns null and the first
    // dereference takes the process down.
    auto params = std::make_unique<SimParam>();
    Param::RegisterParamGroup(params.release());
    sim_ = std::make_unique<Simulation>(TEST_NAME);
  }

  void TearDown() override { sim_.reset(); }

  std::unique_ptr<Simulation> sim_;
};

// case: valid index computes correctly
TEST_F(DiffusionThomasAlgorithmTest, GetBoxIndexValidCoordinates) {
  DiffusionThomasAlgorithm grid(0, "test", 1.0, 0.1, 10, 0.01, false);
  EXPECT_EQ(grid.GetBoxIndex(4, 3, 2), static_cast<size_t>(234));
  // z * res * res + y * res + x = 2*100 + 3*10 + 4 = 234
}

// case: out of bounds triggers assertion
#ifndef NDEBUG
TEST_F(DiffusionThomasAlgorithmTest, GetBoxIndexOutOfBoundsDeath) {
  DiffusionThomasAlgorithm grid(0, "test", 1000.0, 0.1, 10, 0.01, false);
  EXPECT_DEATH(grid.GetBoxIndex(10, 0, 0), "out of bounds");
}
#endif

}  // namespace bdm