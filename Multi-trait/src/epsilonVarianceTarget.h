#ifndef EPSILONVARIANCETARGET_H
#define EPSILONVARIANCETARGET_H
#include <Eigen/Core>
#include "TargetDist.h"

class EpsilonVarianceTarget:public TargetDist{
public:
  //initialize the prior hyper paramters
  EpsilonVarianceTarget(const double _a,const double _b):
  a(_a),
  b(_b){};
  
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
