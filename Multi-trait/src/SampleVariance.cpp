#include <Rcpp.h>
#include <Eigen/Core>
#include "NormalProp.h"
#include "VarianceTarget.h"
// [[Rcpp::depends(RcppEigen)]]
// entry point to sample tha variance from R conditioned on the rest
// [[Rcpp::export]]
Eigen::VectorXd SampleVariance(int idx,const Eigen::MatrixXd& BtB,const Eigen::VectorXd& S, const Eigen::MatrixXd& R, int M,  double v, double tau, int iter, double init ){
  //int _size = R.cols();
  //VarianceTarget target(v, tau,_size );
  //target.R=R;
  //target.S=S;
  //target.BtB=BtB;
  //target.idx = idx;
  //target.M = M;
  //GaussianNewtonian Sampler(&target);
  //Eigen::VectorXd draws_out(iter);
  //Sampler.newtmc_int( init,draws_out,iter,1);
}
