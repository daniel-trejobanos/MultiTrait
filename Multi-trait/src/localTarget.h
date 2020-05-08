// here I target normal(b,tau*c/(tau+c*lambda)) studentT(df,0,1)
#ifndef LOCALTARGET_H
#define LOCALTARGET_H
#include "TargetDist.h"
#include <Eigen/Core>
#include <math.h>

class LocalVarianceTarget : public TargetDist {
public:
  LocalVarianceTarget(double _v0L):
        v0L(_v0L), vp1_2(0.5 * (v0L + 1)),
	lg_vp1_2(std::lgamma(vp1_2)),
        lg_v_2(std::lgamma(0.5 * v0L)),
        logNormTerm(lg_vp1_2 - 0.5 * log(v0L * M_PI) - lg_v_2){};
  void update_log_kernel();
  void update_gradient();
  void update_hessian();
  const double v0L;
  const double vp1_2;
  const double lg_vp1_2;
  const double lg_v_2;
  const double logNormTerm;
  double c;
  double tau;
  double beta2;
};
#endif
