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

Eigen::MatrixXd SamplerGroupVar::sampleGroupVar(int Iter,
                                                const Eigen::VectorXd &BetaSqn,
                                                const Eigen::VectorXd &M0) {

  // group wise squared norms of betas
  target.Bsqnorms = BetaSqn;
  target.m0 = M0;

// if (VerboseTrgt)
#ifdef DEBUG
  printTrgtDiagnostic(target);

  std::cout << "Ngroups" << NGroups << "\n";
  std::cout << "iter" << Iter << "\n";
#endif
  assert(Iter * NGroups > 0);
  Eigen::MatrixXd DrawsOut(Iter, NGroups);
  assert(DrawsOut.size() == Iter * NGroups);
  DrawsOut.setZero();
  Eigen::VectorXd CurrentDraw(10);
  std::cout << "Sampling group variances using Newtonian MC\n";
  for (int i = 0; i < Iter; i++) {
    for (int j = 0; j < NGroups; j++) {
      CurrentDraw.setZero();
      dynamic_cast<GroupVarianceTarget *>(Sampler.target)->idx = j;
      Sampler.newtmc_int(Init(j), CurrentDraw, 1);
      DrawsOut(i, j) = CurrentDraw(9);
    }
    CurrentDraw.setZero();
    DrawsOut(i, NGroups) = CurrentDraw(9);
  }

  return DrawsOut;
}
