#include "NormalProp.h"
#include "dist.h"
#include <Eigen/Core>


// mean of the gaussian distribution
void NormalProp::setParameter1(double point, double gradient, double hessian) {
  parameter1 = point - gradient / hessian;
}

void NormalProp::setParameter2(double point, double gradient, double hessian) {
  parameter2 = -1.0 / hessian;
  assert(parameter2 >= 0);
}

double NormalProp::condLogProb(double conditioned_point) {
  return -0.5 * log(parameter2) -
         0.5 * (conditioned_point - parameter1) / parameter2;
}

double NormalProp::draw() {
  return parameter1 + sqrt(parameter2) * dist.norm_rng(0, 1);
}
