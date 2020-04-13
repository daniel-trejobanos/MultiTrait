#include <iostream>
#include "GammaProp.h"
#include "gammaNewtonian.h"
#include "gaussianNewtonian.h"
#include "groupVarianceTarget.h"
#include <Eigen/Core>
#include "sampleGroupVariances.h"

void printTrgtDiagnostic(GroupVarianceTarget trgt){
  std::cout << "a_sigmaG " << trgt.a << "\n";
  std::cout << "b_sigmaG " << trgt.b << "\n";
  std::cout << "betasqn " << trgt.Bsqnorms << "\n";
  std::cout << "m0 " << trgt.m0 << "\n";
}

Eigen::MatrixXd SamplerGroupVar::sampleGroupVar(int Iter,
                                                 Eigen::VectorXd &BetaSqn,
                                                 Eigen::VectorXd &M0) {
  //target distribution of group sigmaG with gamma priors with parameters A,B
  GroupVarianceTarget target(ASigmaG, BSigmaG);
  //group wise squared norms of betas
  target.Bsqnorms = BetaSqn;
  target.m0 = M0;
  GammaNewtonian Sampler(&target);
  if(VerboseTrgt)
    printTrgtDiagnostic(target);

  Eigen::MatrixXd DrawsOut(Iter, NGroups);
  DrawsOut.setZero();
  Eigen::VectorXd CurrentDraw(10);
  CurrentDraw.setZero();
  
  for (int i = 0; i < Iter; i++) {
    for (int j = 0; j < NGroups; j++) {
      CurrentDraw.setZero();
      dynamic_cast<GroupVarianceTarget*>(Sampler.target)->idx = j;
      Sampler.newtmc_int(Init(j), CurrentDraw, 10, 1);
      DrawsOut(i, j) = CurrentDraw(9);
    }
    CurrentDraw.setZero();
    DrawsOut(i, NGroups) = CurrentDraw(9);
  }

  return DrawsOut;
}
