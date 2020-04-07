#include <Rcpp.h>
#include <Eigen/Core>
#include "GammaProp.h"
#include "gammaNewtonian.h"
#include "gaussianNewtonian.h"
#include "groupVarianceTarget.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// entry point to sample tha variance from R conditioned on the rest
// [[Rcpp::export]]
Eigen::MatrixXd SampleGroupVariance(int iter, 
                               Eigen::VectorXd& betasqn, 
                               Eigen::VectorXd& m0, 
                               double epssqnorm, 
                               double N ,
                               int ngroups,
                               Eigen::VectorXd& a_sigmaG,
                               Eigen::VectorXd& b_sigmaG,
                               double a_eps,
                               double b_eps,
                               double init){
  //
  
  GroupVarianceTarget target(a_sigmaG, b_sigmaG );
  Rcpp::Rcout << "a_sigmaG " << target.a << "\n"; 
  Rcpp::Rcout << "b_sigmaG " << target.b << "\n"; 
  
  target.Bsqnorms = betasqn;
  Rcpp::Rcout << "betasqn " << target.Bsqnorms << "\n"; 
  target.m0 = m0;
  Rcpp::Rcout << "m0 " << target.m0 << "\n"; 
  GammaNewtonian Sampler(&target);
  
    
  Eigen::MatrixXd draws_out(iter,ngroups+1);
  draws_out.setZero();
  Eigen::VectorXd current_draw(10);
  current_draw.setZero();
  for(int i=0; i < iter ; i++){
    for(int j =0; j < ngroups; j++){
      current_draw.setZero();
      target.idx = j;
      Sampler.newtmc_int(0.1,current_draw,10,1);
      draws_out(i,j)= current_draw(9);
    }
    current_draw.setZero();
    //EpsSampler.newtmc_int(R::runif(0,1),current_draw,10,1);
    draws_out(i,ngroups)= current_draw(9);
  }
  
  return draws_out;
}
