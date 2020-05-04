#ifndef SAMPLEPSV_H
#define SAMPLEPSV_H
#include "GammaProp.h"
#include "InvGammaProp.h"
#include "NormalProp.h"
#include "epsilonVarianceTarget.h"
#include "newtonian.h"
#include <Eigen/Core>
class SamplerEpsVar {
public:
  double N;
  const double AEps;
  const double BEps;

  double Init;

  EpsilonVarianceTarget target;

  InvGammaProp prop;
  InvGammaProp prev;

  NewtonianMC Sampler;

  bool VerboseTrgt = 1;
  int seed;
  Distributions_boost &dist;

  SamplerEpsVar(Distributions_boost &dist, double _N, const double _a,
                const double _b, const double _init)
      : dist(dist), N(_N), AEps(_a), BEps(_b), Init(_init), target(_a, _b, _N),
        prop(dist), prev(dist), Sampler(&target, &prop, &prev, dist){};

  double sampleEpsVar(int Iter, const double ESqn);
};
#endif
