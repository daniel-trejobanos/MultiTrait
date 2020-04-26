#include "InvGammaProp.h"
#include "dist.h"
#include "proposal.h"
#include <Eigen/Core>

// mean of the gaussian distribution
void InvGammaProp::setParameter1(double point, double gradient,
                                 double hessian) {
  parameter1 = -2 * ((point * point) * hessian + 2 * point * gradient + 1);

  assert(parameter1 >= 0);
}

void InvGammaProp::setParameter2(double point, double gradient,
                                 double hessian) {
  auto num = abs(point * point * (point * hessian + gradient));

  auto denom = abs((point * point * hessian + 2 * point * gradient + 1));
  auto qut = log(num) - log(denom);

  parameter2 = exp(qut);

  assert(parameter2 >= 0);
}

double InvGammaProp::condLogProb(double conditioned_point) {
  return (0.5 * parameter1) * log(0.5 * parameter1 * parameter2) -
         0.5 * parameter1 * parameter2 / conditioned_point -
         (1 + 0.5 * parameter1) * log(conditioned_point) -
         std::lgamma(0.5 * parameter1);
}

double InvGammaProp::draw() {
  return dist.inv_scaled_chisq_rng(parameter1, parameter2);
}
