#ifndef GAUSSIANNEWTONIAN_H
#define GAUSSIANNEWTONIAN_H
#include <Rcpp.h>
#include "newtonian.h"
#include "NormalProp.h"
#include "TargetDist.h"
using namespace Rcpp;

class GaussianNewtonian: public NewtonianMC{
public:
  GaussianNewtonian(TargetDist *target):NewtonianMC(target){
    NormalProp _prop;
    NormalProp _prev;
    prop = &_prop;
    prev = &_prev; 
  };

};

#endif