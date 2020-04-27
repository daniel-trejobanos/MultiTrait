#ifndef EPSILONVARIANCETARGET_H
#define EPSILONVARIANCETARGET_H
#include "TargetDist.h"
#include <Eigen/Core>

class EpsilonVarianceTarget : public TargetDist {
public:
  // initialize the prior hyper paramters
  EpsilonVarianceTarget(const double _a, const double _b) : a(_a), b(_b){};
  EpsilonVarianceTarget(const double _a, const double _b, double _N)
      : EpsilonVarianceTarget(_a, _b) {
    N = _N;
  };

  void update_log_kernel();
  void update_gradient();
  void update_hessian();

  double Esqnorm;
  double a;
  double b;
  int idx;
  double N;
};
#endif
