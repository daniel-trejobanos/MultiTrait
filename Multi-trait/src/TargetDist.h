#ifndef TARGETDIST_H
#define TARGETDIST_H
#include <Eigen/Core>


class TargetDist{
public:
  TargetDist(){};
  void update(double draw){
    current_value = draw;
    current_value2= pow(current_value,2);
    current_value3=pow(current_value,3);
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
  double current_value2;
  double current_value3;
};
#endif
