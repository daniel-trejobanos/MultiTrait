//
#ifndef VARIANCETARGET_H
#define VARIANCETARGET_H
#include <Eigen/Core>
#include "TargetDist.h"

class VarianceTarget:public TargetDist{
public:
  //initialize the prior hyper paramters
  VarianceTarget(double _v, double _tau,int _size):
  tau(_tau),
  v(_v){
    BtB(_size,_size);
    BtB.setZero();
    S(_size,_size);
    S.setZero();
    R(_size,_size);
    R.setZero();
    Rinv(_size,_size);
    R.setZero();
    idx = 0;
    M = 1;
  };
  
  void update_log_kernel();
  void update_gradient();
  void update_hessian();
  
  Eigen::MatrixXd S;
  Eigen::MatrixXd R;
  Eigen::MatrixXd BtB;
  Eigen::MatrixXd Rinv;
  double v;
  double tau;
  int idx;
  int M;
};
#endif
