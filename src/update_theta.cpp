#include "update_theta.h"
#include <stan/math.hpp>

void update_theta(double &theta, double &acc_theta, const double &mu_theta,
                  const double &jump_theta, const double &sigma,
                  const Eigen::VectorXd &cum_lambda,
                  const Eigen::MatrixXd::ColXpr &beta, const double &gamma,
                  const Eigen::MatrixXd::RowXpr &z, const Eigen::MatrixXd &w,
                  const int &I, const Eigen::MatrixXi::ColXpr &NA,
                  const Eigen::VectorXd &len,
                  const Eigen::MatrixXi::ColXpr &seg,
                  const Eigen::MatrixXd::ColXpr &H,
                  const Eigen::MatrixXi::ColXpr &Y_k, const int cause,
                  boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_theta = 0.0;
  double theta_s;
  // double cum_lambda;

  // std::cout << "Drawing...\n";
  theta_s = stan::math::normal_rng(theta, jump_theta, rng);

  // prior and jump rule
  logr_theta += stan::math::normal_lpdf(theta_s, mu_theta, sigma) -
                stan::math::normal_lpdf(theta, mu_theta, sigma);

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int i = 0; i < I; i++) {

    if (NA(i) == 1) {

      // cum_lambda = 0.0;
      //   for (int g = 0; g < seg(i); g++) {
      //     cum_lambda += len(g) * lambda(i,g);
      //   }
      //   cum_lambda += H(i) * lambda(i,seg(i));

      logr_theta -=
          cum_lambda(i) *
          (stan::math::exp(beta(i) + theta_s -
                           gamma * stan::math::distance(z, w.row(i))) -
           stan::math::exp(beta(i) + theta -
                           gamma * stan::math::distance(z, w.row(i))));

      if (Y_k(i) == cause) {
        logr_theta += theta_s - theta;
      }
    }
  }

  // accept or reject?
  if ((logr_theta > 0.0) ||
      (logr_theta > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    theta = theta_s;
    acc_theta += 1.0;
  }
}
