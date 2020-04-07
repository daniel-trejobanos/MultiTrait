#ifndef TARGETDIST_H
#define TARGETDIST_H
#include <Rcpp.h>
#include <Eigen/Core>


class TargetDist{
public:
  TargetDist(){};
  void update(double draw){current_value = draw;
    update_log_kernel();
    update_gradient();
    update_hessian();};
  virtual void update_log_kernel() =0;
  virtual void update_gradient()=0;
  virtual void update_hessian()=0;
  double log_kernel;
  double gradient;
  double hessian;
  double current_value;
};
#endif
