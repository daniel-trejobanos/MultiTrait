#include "newtonian.h"
#include <Eigen/Core>
#include <random>

double NewtonianMC::newtmc_int(double initial_vals, Eigen::VectorXd &draws_out,
                               int burnin) {

  int iter = draws_out.size();
  bool success = false;

  // setup
  assert(draws_out.size() > 0);
  double first_draw = initial_vals;
#ifdef DEBUG
  std::cout << "first_draw:" << first_draw << "\n";
#endif
  target->update(first_draw);

  double prev_LP = target->log_kernel;
  double prev_gradient = target->gradient;
  double prev_hessian = target->hessian;
  double prop_LP = prev_LP;
  double prop_gradient = target->gradient;
  double prop_hessian = target->hessian;
  double prev_draw = first_draw;
  double new_draw = first_draw;

  // // prev contains the parameers of the proposal in the previous sample
  prev->setParameter1(prev_draw, prev_gradient, prev_hessian);
  prev->setParameter2(prev_draw, prev_gradient, prev_hessian);

#ifdef DEBUG
  std::cout << "prev_LP:" << prev_LP << "\n";
  std::cout << "prev_gradient:" << prev_gradient << "\n";
  std::cout << "prev_hessian:" << prev_hessian << "\n";
  std::cout << "prop_LP:" << prop_LP << "\n";
  std::cout << "prop_gradient:" << prop_gradient << "\n";
  std::cout << "prop_hessian:" << prop_hessian << "\n";
  std::cout << "prev_draw:" << prev_draw << "\n";
  std::cout << "new_draw:" << new_draw << "\n\n";
#endif

  int n_accept = 0;
  double comp_val = 1;
  for (size_t jj = 0; jj < iter; jj++) {
    // we draw a proposal value from the previous parameters
    prop->setParameter1(prev_draw, prev_gradient, prev_hessian);
    prop->setParameter2(prev_draw, prev_gradient, prev_hessian);

#ifdef DEBUG
    std::cout << "\n------------------------------------\n";
    std::cout << "parameter1:" << prop->parameter1 << "\n"
              << "parameter2:" << prop->parameter2 << "\n\n";
#endif

    new_draw = prop->draw();
    // update log prob, gradient and hessian
    target->update(new_draw);
    prop_LP = target->log_kernel;
    prop_gradient = target->gradient;
    prop_hessian = target->hessian;

#ifdef DEBUG
    std::cout << "new_draw:" << new_draw << "\n";
    std::cout << "prop_LP:" << prop_LP << "\n";
    std::cout << "prop_gradient:" << prop_gradient << "\n";
    std::cout << "prop_hessian:" << prop_hessian << "\n\n";
#endif

    // compute acceptance ratio
    comp_val = prop_LP - prev_LP + prop->condLogProb(prev_draw) -
               prev->condLogProb(new_draw);
    double z = dist.unif_rng();

#ifdef DEBUG
    std::cout << "comp_val:" << comp_val << "\n";
    std::cout << "logZ:" << log(z) << "\n\n";
#endif
    //   // Sample if log acc.ratio is bigger than log prob of z
    if (log(z) <= (comp_val)) {
      // now new draw becomes previous draw
      prev_draw = new_draw;
      // we store the proposal log prob, previous gradient and hessian
      prev_LP = prop_LP;
      prev_gradient = prop_gradient;
      prev_hessian = prop_hessian;
      // we also store the parameters of the proposal to compute the
      // conditionals
      prev->parameter1 = prop->parameter1;
      prev->parameter2 = prop->parameter2;
#ifdef DEBUG
      std::cout << "Accepted"
                << "\n\n";
#endif

      draws_out(jj) = new_draw;
      n_accept++;

    } else {

      draws_out(jj) = prev_draw;
    }
  }

  success = true;
  double success_rate = n_accept / (double)iter;
#ifdef DEBUG
  std::cout << "NewtonianMC Acc. Rate:" << success_rate << "\n\n";
#endif
  assert(success_rate > 0);
  return success_rate;
}
