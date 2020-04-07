#ifndef GROUPVARIANCETARGET_H
#define GROUPVARIANCETARGET_H
#include <Rcpp.h>
#include <Eigen/Core>
#include "TargetDist.h"

class GroupVarianceTarget:public TargetDist{
public:
  //initialize the prior hyper paramters
  GroupVarianceTarget(Eigen::VectorXd& _a, Eigen::VectorXd& _b):
  a(_a),
  b(_b){};
  
  void update_log_kernel();
  void update_gradient();
  void update_hessian();
  
  
  Eigen::VectorXd Bsqnorms;
  Eigen::VectorXd m0;
  Eigen::VectorXd sigmaG;
  Eigen::VectorXd a;
  Eigen::VectorXd b;
  int idx;
};
#endif