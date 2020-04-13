#ifndef SAMPLEGV_H
#define SAMPLEGV_H
#include <Eigen/Core>


class SamplerGroupVar {
public:
  double N;
  int NGroups;
  Eigen::VectorXd &ASigmaG;
  Eigen::VectorXd &BSigmaG;
  Eigen::VectorXd &Init;
  bool VerboseTrgt = 0;
  SamplerGroupVar(double _N, int _ngroups, Eigen::VectorXd &_a,
                  Eigen::VectorXd &_b, Eigen::VectorXd& _init)
      : N(_N), NGroups(_ngroups), ASigmaG(_a), BSigmaG(_b), Init(_init){};
  Eigen::MatrixXd sampleGroupVar(int Iter, Eigen::VectorXd &BetaSqn,
                                 Eigen::VectorXd &M0);
};
#endif
