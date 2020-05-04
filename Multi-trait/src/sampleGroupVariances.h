#ifndef SAMPLEGV_H
#define SAMPLEGV_H
#include "InvGammaProp.h"
#include "NormalProp.h"
#include "groupVarianceTarget.h"
#include "newtonian.h"
#include <Eigen/Core>
class SamplerGroupVar {
public:
  double N;
  int NGroups;
  const Eigen::VectorXd &ASigmaG;
  const Eigen::VectorXd &BSigmaG;

  Eigen::VectorXd &Init;

  GroupVarianceTarget target;

  InvGammaProp prop;
  InvGammaProp prev;

  NewtonianMC Sampler;

  bool VerboseTrgt = 0;
  int seed;
  Distributions_boost &dist;

  SamplerGroupVar(Distributions_boost &dist, double _N, int _ngroups,
                  const Eigen::VectorXd &_a, const Eigen::VectorXd &_b,
                   Eigen::VectorXd &_init)
      : dist(dist), N(_N), NGroups(_ngroups), ASigmaG(_a), BSigmaG(_b),
        Init(_init), target(_a, _b), prop(dist), prev(dist),
        Sampler(&target, &prop, &prev, dist){};

  double sampleGroupVar(int Iter, double BetaSqn, double M0, int idx);
};
#endif
