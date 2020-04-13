#include <Rcpp.h>
#include<RcppEigen.h>
#include "sampleGroupVariances.h"
using namespace Rcpp;



// [[Rcpp::export]]
Eigen::MatrixXd GroupSigma(int iter, int N, int ngroups, Eigen::VectorXd& A, Eigen::VectorXd& B, Eigen::VectorXd& init, Eigen::VectorXd& bsqn, Eigen::VectorXd m0 ) {
  SamplerGroupVar Sampler(static_cast<double>(N), ngroups, A, B, init);
  return Sampler.sampleGroupVar(iter,bsqn,m0);
}
