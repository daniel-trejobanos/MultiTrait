#ifndef INVGAMMAPROPOSAL_H
#define INVGAMMAPROPOSAL_H
#include "proposal.h"

class InvGammaProp : public ProposalDist {
public:
  
InvGammaProp(Distributions_boost &_dist) : ProposalDist(_dist) {
  parameter1 = 1; // zero mean
  parameter2 = 1;
};
  void setParameter1(double point, double gradient, double hessian);
  void setParameter2(double pint, double gradient, double hessian);
  double condLogProb(double conditioned_point);
  double draw();
};

#endif
