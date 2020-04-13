//Single site newtonian montecarlo base class of the strategy pattern
#ifndef NEWTONIAN_H
#define NEWTONIAN_H
#include  <Eigen/Core>
#include "proposal.h"
#include "TargetDist.h"

class NewtonianMC{
public:
  NewtonianMC(TargetDist *target):target(target){};
  double newtmc_int(const double initial_vals, Eigen::VectorXd& draws_out, int iter, int burnin);
  ProposalDist *prev; //state of the previous proposal distribution 
  ProposalDist *prop; //state of the current proposal distribution
  TargetDist *target; //state of target distribution
};
#endif
