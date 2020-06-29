#include <stan/math.hpp>
#include "update_theta.h"

void update_theta(double &theta, double &acc_theta, const double &sigma,
                  const Eigen::MatrixXd &lambda,
                  const Eigen::MatrixXd::ColXpr &beta,
                  const double &gamma,
                  const Eigen::MatrixXd::RowXpr &z, const Eigen::MatrixXd &w,
                  const int &I,
                  const Eigen::VectorXd &len, const Eigen::MatrixXi::ColXpr &seg,
                  const Eigen::MatrixXd::ColXpr &H,
                  const Eigen::MatrixXi::ColXpr &Y_k, const int cause,
                  boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_theta = 0.0;
  double theta_s;
  double cum_lambda;

  // std::cout << "Drawing...\n";
  theta_s = stan::math::normal_rng(theta, 1.0, rng);

  // prior and jump rule
  logr_theta +=
      stan::math::normal_lpdf(theta_s, 0.0, sigma) -
      stan::math::normal_lpdf(theta, 0.0, sigma);

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int i = 0; i < I; i++) {
  cum_lambda = 0.0;
    for (int g = 0; g < seg(i); g++) {
      cum_lambda += len(g) * lambda(i,g);
    }
    cum_lambda += H(i) * lambda(i,seg(i));

    logr_theta -= cum_lambda * ( stan::math::exp( beta(i) + theta_s - gamma * stan::math::distance(z, w.row(i))) -
                                 stan::math::exp( beta(i) + theta - gamma * stan::math::distance(z, w.row(i))));

      if (Y_k(i) == cause) {
        logr_theta += stan::math::log(theta_s) - stan::math::log(theta);
      }
    }

  // accept or reject?
  if ((logr_theta) > 0.0 || (logr_theta >
                              stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    theta = theta_s;
    acc_theta += 1.0;
  }
}


