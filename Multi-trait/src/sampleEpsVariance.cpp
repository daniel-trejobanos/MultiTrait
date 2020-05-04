#include "sampleEpsVariance.h"
#include "epsilonVarianceTarget.h"
#include <Eigen/Core>
#include <iostream>

void printTrgtDiagnostic(EpsilonVarianceTarget trgt) {
  std::cout << "a_eps: " << trgt.a << "\n\n";
  std::cout << "b_eps: " << trgt.b << "\n\n";
  std::cout << "Esqn: " << trgt.Esqnorm << "\n\n";
  std::cout << "N " << trgt.N << "\n\n";
}

double SamplerEpsVar::sampleEpsVar(int Iter,const double ESqn) {
  // target distribution of group sigmaG with gamma priors with parameters A,B
  // GroupVarianceTarget target(ASigmaG, BSigmaG);
  // group wise squared norms of betas
  target.Esqnorm = ESqn;
   target.N = N;
  // if (VerboseTrgt)
  printTrgtDiagnostic(target);
  std::cout<< "iter" <<Iter<<"\n";
  Eigen::VectorXd DrawsOut(Iter);
  DrawsOut.setZero();
  Eigen::VectorXd CurrentDraw(10);
  //Init=dist.unif_rng();
  for (int i = 0; i < Iter; i++) {
    
    CurrentDraw.setZero();
    Sampler.newtmc_int(Init, CurrentDraw, 1);
    
  }
  //to ease integration we only return the last iteration
  return CurrentDraw(9);
}
