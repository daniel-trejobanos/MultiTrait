#include <Eigen/Core>
#include "groupVarianceTarget.h"

double inline priorlpdf(double x, double a, double b){
  return a*log(b) + (a -1 ) * log(x) - b*x - std::lgamma(a);
}

double inline gradientPrior(double x, double a, double b){
  return  -b + ( a-1 ) / x;
}

double inline hessianPrior(double x, double a, double b){
  return -(a-1)/(x*x);
}


void GroupVarianceTarget::update_log_kernel(){
  double log_lik =  -0.5 * m0(idx) * Bsqnorms(idx) / current_value - 0.5*log(current_value);
  log_kernel = log_lik + priorlpdf(current_value, a(idx), b(idx));
}


void GroupVarianceTarget::update_gradient(){
  auto gradient_log_lik = 0.5 * m0(idx) * Bsqnorms(idx) / (current_value2) - 0.5/current_value;
  gradient = gradient_log_lik + gradientPrior(current_value, a(idx), b(idx));
}


void GroupVarianceTarget::update_hessian(){
  auto hessian_log_lik = - m0(idx) * Bsqnorms(idx) / (current_value3) + 0.5/(current_value2);
  hessian= hessian_log_lik + hessianPrior(current_value,a(idx),b(idx));
}
