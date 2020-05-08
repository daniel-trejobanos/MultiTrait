#include "localTarget.h"
#include <Eigen/Core>

void LocalVarianceTarget::update_log_kernel() {
  auto priorlpdf = logNormTerm - vp1_2 * log(1 + current_value2 / v0L);
  auto lambda_tilde = c * tau / (tau + c * current_value);
  auto log_lik = -0.5 * log(lambda_tilde) - 0.5 * beta2 / lambda_tilde;
  log_kernel = log_lik + priorlpdf;
}

void LocalVarianceTarget::update_gradient() {
  auto gradientPrior =
      -vp1_2 * 2 * current_value / (v0L * (1 + current_value2 / v0L));
  auto gradient_log_lik = 0.5 * (-beta2 / tau + c / (tau + c * current_value));

  gradient = gradient_log_lik + gradientPrior;
}

void LocalVarianceTarget::update_hessian() {
  auto nat = 1 + current_value2 / v0L;
  auto hessianPrior =
      -vp1_2 * (4 * current_value2 / (v0L * v0L * nat * nat) + 2 / (v0L * nat));
  auto hessian_log_lik =
      -0.5 * c * c / ((tau + c * current_value) * (tau + c * current_value));
  hessian = hessian_log_lik + hessianPrior;
}

