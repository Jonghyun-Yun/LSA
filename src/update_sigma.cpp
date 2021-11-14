#include "update_sigma.h"

double update_sigma(const double& a_sigma, const double &b_sigma,
                           const Eigen::MatrixXd &theta, const int &N, boost::ecuyer1988 &rng) {
  return stan::math::sqrt(stan::math::inv_gamma_rng(a_sigma + 0.5 * theta.array().size(), b_sigma + 0.5 * theta.array().square().sum(), rng));
}
