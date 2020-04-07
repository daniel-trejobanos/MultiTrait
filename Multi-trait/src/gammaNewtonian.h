#ifndef GAMMANEWTONIAN_H
#define GAMMANEWTONIAN_H
#include <Rcpp.h>
#include "newtonian.h"
#include "GammaProp.h"
#include "TargetDist.h"
using namespace Rcpp;

class GammaNewtonian: public NewtonianMC{
public:
  GammaNewtonian(TargetDist *target):NewtonianMC(target){
    GammaProp _prop;
    GammaProp _prev;
    prop = &_prop;
    prev = &_prev; 
  };
  
};

#endif
