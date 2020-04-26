#ifndef NORMALPROPOSAL_H
#define NORMALPROPOSAL_H
#include "proposal.h"

class NormalProp : public ProposalDist {
public:
  NormalProp(Distributions_boost &dist):ProposalDist(dist){parameter1=0;parameter2=1;};
  void setParameter1(double point, double gradient, double hessian);
  void setParameter2(double pint, double gradient, double hessian);
  double condLogProb(double conditioned_point);
  double draw();
};

#endif
