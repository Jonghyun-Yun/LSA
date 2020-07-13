#include <stan/math.hpp>
#include "update_beta.h"

Eigen::VectorXd update_beta(double &beta, double &acc_beta,
                           const double &mu_beta, const double &sigma_beta, const double &jump_beta,
                           const Eigen::MatrixXd::RowXpr &lambda,
                           const Eigen::MatrixXd::ColXpr &theta,
                           const double &gamma,
                           const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
                           const int &N,
                           const Eigen::MatrixXi::RowXpr &NA,
                           const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                           const Eigen::MatrixXd::RowXpr &H,
                           const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                           boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_beta = 0.0;
  double beta_s;
  Eigen::VectorXd cum_lambda = Eigen::VectorXd::Zero(N);

  // // std::cout << "Drawing...\n";
  beta_s = stan::math::normal_rng(beta, jump_beta, rng);

  // prior and jump rule
  logr_beta +=
      stan::math::normal_lpdf(beta_s, mu_beta, sigma_beta) -
      stan::math::normal_lpdf(beta, mu_beta, sigma_beta);

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int k = 0; k < N; k++) {

    if (NA(k) == 1) {

      for (int g = 0; g < seg(k); g++) {
        cum_lambda(k) += len(g) * lambda(g);
      }
      cum_lambda(k) += H(k) * lambda(seg(k));

      logr_beta -= cum_lambda(k) * ( stan::math::exp( beta_s + theta(k) - gamma * stan::math::distance(z.row(k), w)) -
                                     stan::math::exp( beta + theta(k) - gamma * stan::math::distance(z.row(k), w)));

      if (Y_i(k) == cause) {
        logr_beta += stan::math::log(beta_s) - stan::math::log(beta);
      }

    }
  }
 
  // accept or reject?
  if ((logr_beta > 0.0) ||
      (logr_beta > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    beta = beta_s;
    acc_beta += 1.0;
  }

  return cum_lambda;
}


