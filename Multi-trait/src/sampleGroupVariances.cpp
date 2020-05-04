#include "sampleGroupVariances.h"
#include "groupVarianceTarget.h"
#include <Eigen/Core>
#include <iostream>

void printTrgtDiagnostic(GroupVarianceTarget trgt) {
  std::cout << "a_sigmaG " << trgt.a << "\n";
  std::cout << "b_sigmaG " << trgt.b << "\n";
  std::cout << "betasqn " << trgt.Bsqnorms << "\n";
  std::cout << "m0 " << trgt.m0 << "\n";
}

double SamplerGroupVar::sampleGroupVar(int Iter, double BetaSqn, double M0, int idx) {

  // group wise squared norms of betas
  dynamic_cast<GroupVarianceTarget *>(Sampler.target)->Bsqnorms = BetaSqn;
  dynamic_cast<GroupVarianceTarget *>(Sampler.target)->m0= M0;
  dynamic_cast<GroupVarianceTarget *>(Sampler.target)->idx = idx;
  // newtonian seems very sensitive to initial conditions
  // Init(idx) = dist.unif_rng();

  // if (VerboseTrgt)

  printTrgtDiagnostic(target);

  std::cout << "Ngroups" << NGroups << "\n";
  std::cout << "iter" << Iter << "\n";

  assert(Iter * NGroups > 0);
  Eigen::MatrixXd DrawsOut(Iter, NGroups);
  assert(DrawsOut.size() == Iter * NGroups);
  DrawsOut.setZero();
  Eigen::VectorXd CurrentDraw(10);
  std::cout << "Sampling group variances using Newtonian MC\n";
  for (int i = 0; i < Iter; i++) {

    CurrentDraw.setZero();

    Sampler.newtmc_int(Init(idx), CurrentDraw, 1);
  }


return CurrentDraw(9);
}
