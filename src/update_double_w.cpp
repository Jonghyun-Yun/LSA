#include "update_double_w.h"
#include <stan/math.hpp>

Eigen::VectorXd update_double_w(const Eigen::MatrixXd::RowXpr &w, double &acc_w,
                                const Eigen::MatrixXd::RowXpr &mu_w,
                                const Eigen::MatrixXd::RowXpr &sigma_w,
                                const Eigen::MatrixXd::RowXpr &jump_w,
                                const Eigen::VectorXd &cum_lambda,
                                const double &beta,
                                const Eigen::MatrixXd::ColXpr &theta,
                                const double &gamma,
                                const Eigen::MatrixXd &z,
                                const int &N,
                                const Eigen::MatrixXi::RowXpr &NA,
                                const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                                const Eigen::MatrixXd::RowXpr &H,
                                const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                                boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_w = 0.0;
  Eigen::VectorXd w_s(2);

  // std::cout << "Drawing...\n";
  for (int d = 0; d < 2; d++)
  {
      w_s(d) = stan::math::normal_rng(w(d), jump_w(d), rng);
      logr_w += stan::math::normal_lpdf(w_s(d), mu_w(d), sigma_w(d)) -
          stan::math::normal_lpdf(w(d), mu_w(d), sigma_w(d));
  }

    // w_s = stan::math::normal_rng(w, jump_w, rng);
    // logr_w +=
    //   stan::math::normal_lpdf(w_s, mu_w, sigma_w) -
    //   stan::math::normal_lpdf(w, mu_w, sigma_w);

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int k = 0; k < N; k++) {

    if (NA(k) == 1) {

        logr_w -=
          cum_lambda(k) *
          (stan::math::exp(beta + theta(k) -
                           gamma * stan::math::distance(z.row(k), w_s)) -
           stan::math::exp(beta + theta(k) -
                           gamma * stan::math::distance(z.row(k), w)));

        if (Y_i(k) == cause) {
          logr_w -= gamma * (stan::math::distance(z.row(k), w_s) -
                                stan::math::distance(z.row(k), w));
        }
      }
    }

  // accept or reject?
  if ((logr_w > 0.0) ||
      (logr_w > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    acc_w += 1.0;
  } else {
    w_s = w; // do NOT aceept w_s, use the current w instead
  }
  return w_s;
}
