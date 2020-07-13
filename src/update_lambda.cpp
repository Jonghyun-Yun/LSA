#include <stan/math.hpp>
#include "update_lambda.h"

void update_lambda(double &lambda, double &acc_lambda, const double &a_lambda,
                   const double &b_lambda, const double &jump_lambda,
                   const int g, const double &beta,
                   const Eigen::MatrixXd::ColXpr &theta, const double &gamma,
                   const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
                   const int &N, const Eigen::MatrixXi::RowXpr &NA,
                   const double &len, const Eigen::MatrixXi::RowXpr &seg,
                   const Eigen::MatrixXd::RowXpr &H,
                   const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                   boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_lambda = 0.0;
  double lambda_s;

  if (NA.array().sum() == 0) {

    lambda = stan::math::gamma_rng(a_lambda, b_lambda, rng);
    acc_lambda += 1.0;

  } else {

    // std::cout << "Drawing...\n";
    lambda_s =
        stan::math::lognormal_rng(stan::math::log(lambda), jump_lambda, rng);

    // prior and jump rule
    logr_lambda += stan::math::gamma_lpdf(lambda_s, a_lambda, b_lambda) -
                   stan::math::gamma_lpdf(lambda, a_lambda, b_lambda) +
                   stan::math::lognormal_lpdf(lambda, stan::math::log(lambda_s),
                                              jump_lambda) -
                   stan::math::lognormal_lpdf(lambda_s, stan::math::log(lambda),
                                              jump_lambda);

    // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
    for (int k = 0; k < N; k++) {

      if (NA(k) == 1) {

        if (g < seg(k)) {
          logr_lambda -=
              len * (lambda_s - lambda) *
              stan::math::exp(beta + theta(k) -
                              gamma * stan::math::distance(z.row(k), w));
        } else if (g == seg(k)) {
          logr_lambda -=
              H(k) * (lambda_s - lambda) *
              stan::math::exp(beta + theta(k) -
                              gamma * stan::math::distance(z.row(k), w));
          if (Y_i(k) == cause) {
            logr_lambda += stan::math::log(lambda_s) - stan::math::log(lambda);
          }
        }
      }
    }

    // accept or reject?
    if ((logr_lambda > 0.0) ||
        (logr_lambda >
         stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
      lambda = lambda_s;
      acc_lambda += 1.0;
    }
  }
}
