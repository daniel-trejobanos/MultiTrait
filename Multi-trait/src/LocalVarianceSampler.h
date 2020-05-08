#ifndef SAMPLERLOCAL_H
#define SAMPLERLOCAL_H
#include "InvGammaProp.h"
#include "localTarget.h"
#include "newtonian.h"
#include <Eigen/Core>

class SamplerLocalVar {

public:
  LocalVarianceTarget target;
  InvGammaProp prop;
  InvGammaProp prev;
  NewtonianMC Sampler;
  double Init;
  const double v0L;

  Distributions_boost &dist;

  SamplerLocalVar(Distributions_boost &dist, const double _v0L,
                  const double _init)
    : dist(dist), Init(_init), v0L(_v0L), target(v0L),prop(dist),prev(dist),Sampler(&target, &prop, &prev, dist){};
  double sampleLocalVar(int Iter, double tau, double c, double beta2);
};

#endif
