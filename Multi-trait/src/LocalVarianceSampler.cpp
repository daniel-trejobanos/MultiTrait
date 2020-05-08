// TODO
#include "LocalVarianceSampler.h"

double SamplerLocalVar::sampleLocalVar(int Iter, double tau, double c,
                                       double beta2) {
  
  dynamic_cast<LocalVarianceTarget *>(Sampler.target)->beta2 = beta2;
  dynamic_cast<LocalVarianceTarget *>(Sampler.target)->c = c;
  dynamic_cast<LocalVarianceTarget *>(Sampler.target)->tau = tau;
  Eigen::VectorXd DrawsOut(Iter);
  DrawsOut.setZero();
  Sampler.newtmc_int(Init, DrawsOut, 1);

  return DrawsOut(9);
}
