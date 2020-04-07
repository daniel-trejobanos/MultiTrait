#ifndef GAMMAPROPOSAL_H
#define GAMMAPROPOSAL_H
#include "proposal.h"

class GammaProp:public ProposalDist{
public:
  GammaProp();
  void setParameter1(double point, double gradient, double hessian);
  void setParameter2(double pint, double gradient, double hessian);
  double condLogProb(double conditioned_point);
  double draw();
};

#endif