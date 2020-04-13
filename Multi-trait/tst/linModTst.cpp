#include "GammaProp.h"
#include "gammaNewtonian.h"
#include "groupVarianceTarget.h"
#include "linearModel.h"
#include "sampleGroupVariances.h"
#include "gtest/gtest.h"
#include <math.h>

namespace my {
namespace project {
namespace {
class LinearModelSimulator : public ::testing::Test {
protected:
  void SetUp() override {
    LM = LinearModel(4000, 100, 2000, 3, 0.6);
    LM.simulate();
  }
  LinearModel LM;
};

typedef LinearModelSimulator LinearModelTester;
typedef LinearModelSimulator GVSamplerTester;

TEST_F(LinearModelTester, phenMean) {
  double mean = LM.Y.mean();
  EXPECT_DOUBLE_EQ(mean, 0);
}

TEST_F(LinearModelTester, phenSd) {
  double var = (Eigen::VectorXd(LM.Y.array() - LM.Y.mean())).squaredNorm() /
               (LM.Y.size());
  EXPECT_DOUBLE_EQ(var, 1);
}

TEST_F(LinearModelTester, EMean) {
  double mean = LM.E.mean();
  EXPECT_DOUBLE_EQ(mean, 0);
}

TEST_F(LinearModelTester, ESd) {
  double var = (Eigen::VectorXd(LM.E.array() - LM.E.mean())).squaredNorm() /
               (LM.E.size());
  EXPECT_DOUBLE_EQ(var, 1.0 - LM.h2);
}

TEST_F(GVSamplerTester, testLogKernel) {

  double N = static_cast<double>(LM.Y.rows());
  int Ngroups = LM.NumGroups;
  Eigen::VectorXd A(Ngroups);    // a parameter of prior
  Eigen::VectorXd B(Ngroups);    // b parameter of prior
  double Init = LM.h2 / Ngroups; // start at true value
  Eigen::VectorXd M0(Ngroups);
  M0.array() = ceil(LM.M / Ngroups);
  A.array() = 2 / 3;
  B.array() = 2;
  GroupVarianceTarget target(A, B);
  // group wise squared norms of betas
  target.Bsqnorms = LM.BSqN;
  target.m0 = M0;
  target.idx = 0;
  target.update(Init);
  // TODO for loop to test all log kernels
  double llk = -0.5 * M0(0) * LM.BSqN(0) / Init - 0.5 * log(Init) -
               A(0) * log(B(0)) + (A(0) - 1) * log(Init) - B(0) * Init -
               std::lgamma(A(0));
  ASSERT_DOUBLE_EQ(llk, target.log_kernel);
}

TEST_F(GVSamplerTester, Testsamplerlogkernel) {

  double N = static_cast<double>(LM.Y.rows());
  int Ngroups = LM.NumGroups;
  Eigen::VectorXd A(Ngroups); // a parameter of prior
  Eigen::VectorXd B(Ngroups); // b parameter of prior
  Eigen::VectorXd Init(Ngroups);
  Init.array() = LM.h2 / Ngroups; // start at true value
  Eigen::VectorXd M0(Ngroups);
  M0.array() = ceil(LM.M / Ngroups);
  A.array() = 2 / 3;
  B.array() = 2;
  GroupVarianceTarget target(A, B);

  // group wise squared norms of betas
  target.Bsqnorms = LM.BSqN;
  target.m0 = M0;
  GammaNewtonian Sampler(&target);
  dynamic_cast<GroupVarianceTarget *>(Sampler.target)->idx = 0;
  Sampler.target->update(Init(0));
  // TODO for loop to test all log kernels
  double llk = -0.5 * M0(0) * LM.BSqN(0) / Init(0) - 0.5 * log(Init(0)) -
               A(0) * log(B(0)) + (A(0) - 1) * log(Init(0)) - B(0) * Init(0) -
               std::lgamma(A(0));
  ASSERT_DOUBLE_EQ(llk, Sampler.target->log_kernel);
}


} // Namespace
} // namespace project
} // namespace my
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
