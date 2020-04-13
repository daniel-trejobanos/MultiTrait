#include <Eigen/Core>
#include "GammaProp.h"
#include "dist.h"


GammaProp::GammaProp(){
  parameter1 = 1; // zero mean
  parameter2 = 1;
  dist= Distributions_boost(123);
}

//mean of the gaussian distribution
void GammaProp::setParameter1(double point, double gradient, double hessian){
  parameter1 = 1 -point* point * hessian;
  assert(parameter1 >= 0 );
}

void GammaProp::setParameter2(double point, double gradient, double hessian){
  parameter2 = - point *hessian -gradient;
  assert(parameter2 >= 0);
}

double GammaProp::condLogProb(double conditioned_point){
  return parameter1*log(parameter2) + (parameter1 -1 ) * log(conditioned_point) - parameter2*conditioned_point - std::lgamma(parameter1);
}

double GammaProp::draw(){
  return dist.rgamma(parameter1, parameter2) ;
}
