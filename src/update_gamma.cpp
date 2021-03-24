#include <stan/math.hpp>
#include "update_gamma.h"

void update_gamma(Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
                  const Eigen::VectorXd &mu_gamma,
                  const Eigen::VectorXd &sigma_gamma,
                  const Eigen::VectorXd &jump_gamma,
                  const Eigen::MatrixXd &cum_lambda,
                  const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta,
                  const Eigen::MatrixXd &z, const Eigen::MatrixXd &w,
                  const int &I, const int &N, const int &G,
                  const Eigen::MatrixXi &NA,
                  const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                  const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y,
                  boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  Eigen::VectorXd logr_gamma(2); // gamma for each cause
  logr_gamma.setZero();
  Eigen::VectorXd gamma_s(2);
  double gamma_c; // when gamma(c) < 0
  // Eigen::VectorXd cum_lambda(2); // cumulative lambda for each cause
  // Eigen::MatrixXd lambda(I, G);
  // Eigen::MatrixXd z(N, 2);

  for (int c = 0; c < 2; c++) {

    // if (c == 0) {
    //   lambda = lambda0;
    //   z = z0;
    // } else {
    //   lambda = lambda1;
    //   z = z1;
    // }

    if (gamma(c) < 0) {
      gamma_c = abs(gamma(c)); // pretend gamma(c) is positive
    // std::cout << "Drawing...\n";
    gamma_s(c) = stan::math::lognormal_rng(stan::math::log(gamma_c), jump_gamma(c), rng);

    // prior and jumping rule
    logr_gamma(c) +=
      stan::math::lognormal_lpdf(gamma_s(c), mu_gamma(c), sigma_gamma(c)) -
      stan::math::lognormal_lpdf(gamma_c, mu_gamma(c), sigma_gamma(c)) +
      stan::math::lognormal_lpdf(gamma_c, stan::math::log(gamma_s(c)), jump_gamma(c)) -
      stan::math::lognormal_lpdf(gamma_s(c), stan::math::log(gamma_c), jump_gamma(c));
    gamma_s(c) *= -1.0; // put it back to negative
    } else {
    // std::cout << "Drawing...\n";
    gamma_s(c) = stan::math::lognormal_rng(stan::math::log(gamma(c)), jump_gamma(c), rng);

    // prior and jumping rule
    logr_gamma(c) +=
      stan::math::lognormal_lpdf(gamma_s(c), mu_gamma(c), sigma_gamma(c)) -
      stan::math::lognormal_lpdf(gamma(c), mu_gamma(c), sigma_gamma(c)) +
      stan::math::lognormal_lpdf(gamma(c), stan::math::log(gamma_s(c)), jump_gamma(c)) -
      stan::math::lognormal_lpdf(gamma_s(c), stan::math::log(gamma(c)), jump_gamma(c));
    }
    // std::cout << "Calculating the log-acceptance ratio...\n";
    for (int i = 0; i < I; i++) {
      for (int k = 0; k < N; k++) {

        if (NA(i,k) == 1) {

          // cum_lambda.setZero();

          // for (int g = 0; g < seg(i,k); g++) {
          //   cum_lambda(c) += len(g) * lambda(i,g);
          // }
          // cum_lambda(c) += H(i,k) * lambda(i,seg(i,k));

          logr_gamma(c) -=
            cum_lambda(c * I + i, k) *
            (stan::math::exp(beta(i, c) + theta(k, c) -
                             gamma_s(c) *
                             stan::math::distance(z.row(c*N + k), w.row(c*I + i))) -
             stan::math::exp(beta(i, c) + theta(k, c) -
                             gamma(c) *
                             stan::math::distance(z.row(c*N + k), w.row(c*I + i))));

          if (Y(i, k) == c) {
            logr_gamma(c) -= (gamma_s(c) - gamma(c)) * stan::math::distance(z.row(c*N + k), w.row(c*I + i));
          }
        }
      }
    }

    // accept or reject?
    if ((logr_gamma(c) > 0.0) || (logr_gamma(c) >
                                  stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
      acc_gamma(c) += 1.0;
      gamma(c) = gamma_s(c);
    }
  }
}
