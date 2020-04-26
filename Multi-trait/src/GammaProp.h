#ifndef GAMMAPROPOSAL_H
#define GAMMAPROPOSAL_H
#include "proposal.h"

class GammaProp : public ProposalDist {
public:
  
GammaProp(Distributions_boost &_dist) : ProposalDist(_dist) {
  parameter1 = 1; // zero mean
  parameter2 = 1;
};
  void setParameter1(double point, double gradient, double hessian);
  void setParameter2(double pint, double gradient, double hessian);
  double condLogProb(double conditioned_point);
  double draw();
};

#endif
