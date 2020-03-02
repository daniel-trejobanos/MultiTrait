#ifndef NORMALPROPOSAL_H
#define NORMALPROPOSAL_H
#include "proposal.h"
using namespace Rcpp;

class NormalProp:public ProposalDist{
public:
  NormalProp();
  void setParameter1(double point, double gradient, double hessian);
  void setParameter2(double pint, double gradient, double hessian);
  double condLogProb(double conditioned_point);
  double draw();
};

#endif