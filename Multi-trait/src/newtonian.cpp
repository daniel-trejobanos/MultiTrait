#include <Rcpp.h>
#include <RcppEigen.h>
#include "newtonian.h"
#include "NormalProp.h"
using namespace Rcpp;
double
  NewtonianMC::newtmc_int(const double initial_vals, Eigen::VectorXd& draws_out, int iter, int burnin)
  {
   
    bool success = false;
    
    
      // setup
      
     double first_draw = initial_vals;
      
      
      target->update(first_draw);
      
      auto prev_LP = target->log_kernel;
      auto prev_gradient = target->gradient;
      auto prev_hessian = target->hessian;
      auto prop_LP = prev_LP;
      auto prop_gradient = target->gradient;
      auto prop_hessian = target->hessian;
      
      auto prev_draw = first_draw;
      auto new_draw  = first_draw;
    
      //prev contains the parameers of the proposal in the previous sample
      prev->setParameter1(prev_draw, target->gradient, target->hessian);
      prev->setParameter2(prev_draw, target->gradient, target->hessian);
  
      
      //
      
      int n_accept = 0;
      
      for (size_t jj = 0; jj < draws_out.size(); jj++)
      {
        //we draw a proposal value from the previous parameters
        //
        prop->setParameter1(prev_draw, prev_gradient, prev_hessian);
        prop->setParameter2(prev_draw, prev_gradient, prev_hessian);
        new_draw = prop->draw();
        //update log prob, gradient and hessian
        target->update(new_draw);
        prop_LP = target->log_kernel;
        prop_gradient = target->gradient;
        prop_hessian = target->hessian;
        
        //
        //compute acceptance ratio
        double comp_val = std::min(0.0, prop_LP - prev_LP + prop->condLogProb(prev_draw) - prev->condLogProb(new_draw));
        double z = R::runif(0,1);
        //sample if log acc.ratio is bigger than log prob of z
        if (log(z) <= (comp_val))
        {
          // now new draw becomes previous draw
          prev_draw = new_draw;
          //we store the proposal log prob, previous gradient and hessian
          prev_LP = prop_LP;
          prev_gradient = prop_gradient;
          prev_hessian = prop_hessian;
          //we also store the parameters of the proposal to compute the conditionals
          prev->parameter1 = prop->parameter1;
          prev->parameter2 = prop->parameter2;
          
          draws_out(jj) = new_draw;
          n_accept++;
          
        }
        else
        {
          
          draws_out(jj) = prev_draw;
          
        }
      }
      
      success = true;
      auto success_rate = n_accept/draws_out.size();
   
      
      return success_rate;
  }
