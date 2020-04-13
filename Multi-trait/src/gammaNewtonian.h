#ifndef GAMMANEWTONIAN_H
#define GAMMANEWTONIAN_H
#include "newtonian.h"
#include "GammaProp.h"
#include "TargetDist.h"

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
