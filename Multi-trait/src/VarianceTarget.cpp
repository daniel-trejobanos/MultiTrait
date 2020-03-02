#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include "VarianceTarget.h"

double inline priorlpdf(double x, double v, double tau){
  return - v*tau/(2*x) - ( 1 + 0.5*v ) * log(x);
}

double inline gradientPrior(double x, double v, double tau){
  return  v*tau/(2*x*x) - ( 1 + 0.5*v ) / x;
}

double inline hessianPrior(double x, double v, double tau){
  return - 2* v * tau/(2*x*x*x) + ( 1+0.5 * v ) /(x*x);
}

//function to compute the dot product of the offdiagonal elements
double off_diagonal_product(const Eigen::VectorXd& A, const Eigen::VectorXd& B, const Eigen::VectorXd& C, int idx ){
  double off_diagonal_prod=0;
  for(int ii = 0; ii< A.size() ; ii++)
  {
    off_diagonal_prod+= A(ii) * B(ii) / C(ii);
  }
  return off_diagonal_prod - A(idx) * B(idx) /C(idx);
}
//s;ow form of the updates, just while I test.
void VarianceTarget::update_log_kernel(){
  auto log_variance =  - 0.5 * M * log(current_value);
  auto off_diagonal = off_diagonal_product(BtB.row(idx), Rinv.col(idx), S, idx) /sqrt(current_value);
  auto diagonal = BtB(idx,idx) * Rinv(idx,idx) / current_value; 
  auto log_trace =   -off_diagonal -0.5* diagonal;
  log_kernel =  log_variance + log_trace  +  priorlpdf(current_value,v, tau);
}


void VarianceTarget::update_gradient(){
  auto gradient_log_variance = - 0.5 * M / current_value;
  auto gradient_diagonal =  -Rinv(idx,idx) * BtB(idx,idx) /(current_value*current_value);
  auto gradient_offdiagonal = -0.5*off_diagonal_product(BtB.row(idx), Rinv.col(idx), S, idx) / sqrt(current_value * current_value * current_value);
  auto gradient_log_trace= -gradient_offdiagonal - 0.5*gradient_diagonal;
  gradient = gradient_log_variance  + gradient_log_trace + gradientPrior(current_value, v, tau);
}


void VarianceTarget::update_hessian(){
  auto hessian_log_variance =  0.5 * M / (current_value*current_value);
  auto hessian_diagonal =  2 * Rinv(idx,idx) * BtB(idx,idx) /(current_value*current_value * current_value);
  auto hessian_offdiagonal = (1.5)*off_diagonal_product(BtB.row(idx), Rinv.col(idx), S, idx) / sqrt(pow(current_value,  5));
  auto hessian_log_trace= -hessian_offdiagonal - 0.5*hessian_diagonal;
  hessian =  hessian_log_variance + hessian_log_trace + hessianPrior(current_value, v ,tau);
}