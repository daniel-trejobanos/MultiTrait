// proposal distribution abstract class, each proposal must have two parameters, the 
#ifndef PROPOSAL_H
#define PROPOSAL_H
#include "dist.h"

class ProposalDist{
public:
  //constructor with a given proposal at a given point, no draws possible before initialising into a point.
  virtual void setParameter1(double point, double gradient, double hessian) = 0;
  virtual void setParameter2(double pint, double gradient, double hessian) = 0;
  virtual double condLogProb(double conditioned_point) =0;
  virtual double draw()=0;
 
  double parameter1;
  double parameter2;
  Distributions_boost dist;
};

#endif
