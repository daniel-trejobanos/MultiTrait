// Single site newtonian montecarlo base class of the strategy pattern
#ifndef NEWTONIAN_H
#define NEWTONIAN_H
#include "TargetDist.h"
#include "proposal.h"
#include <Eigen/Core>

class NewtonianMC {
public:
  NewtonianMC(TargetDist *target, ProposalDist *prop, ProposalDist *prev,
              Distributions_boost &dist)
      : target(target), prev(prev), prop(prop), dist(dist){};
  double newtmc_int(const double initial_vals, Eigen::VectorXd &draws_out,
                    int burnin);
  ProposalDist *prev; // state of the previous proposal distribution
  ProposalDist *prop; // state of the current proposal distribution
  TargetDist *target; // state of target distribution
  const bool showDebug = 1;
  Distributions_boost &dist;
};
#endif
