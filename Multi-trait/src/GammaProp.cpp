#include <Rcpp.h>
#include <Eigen/Core>
#include "GammaProp.h"
using namespace Rcpp;

GammaProp::GammaProp(){
  parameter1 = 1; // zero mean
  parameter2 = 1;
}

//mean of the gaussian distribution
void GammaProp::setParameter1(double point, double gradient, double hessian){
  parameter1 = 1 -point* point * hessian;
}

void GammaProp::setParameter2(double point, double gradient, double hessian){
  parameter2 = - point *hessian -gradient;
}

double GammaProp::condLogProb(double conditioned_point){
  return parameter1*log(parameter2) + (parameter1 -1 ) * log(conditioned_point) - parameter2*conditioned_point - std::lgamma(parameter1);
}

double GammaProp::draw(){
  return R::rgamma(parameter1, 1/parameter2) ;
}
