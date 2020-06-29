#include <stan/math.hpp>
#include "update_w.h"

Eigen::VectorXd update_w(const Eigen::MatrixXd::RowXpr &w, double &acc_w,
                         const Eigen::MatrixXd::RowXpr &lambda0,
                         const Eigen::MatrixXd::RowXpr &lambda1,
                         const Eigen::MatrixXd::RowXpr &beta,
                         const Eigen::MatrixXd &theta,
                         const Eigen::VectorXd &gamma,
                         const Eigen::MatrixXd &z0, const Eigen::MatrixXd &z1,
                         const int &N, const int &G, const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                         const Eigen::MatrixXd::RowXpr &H, const Eigen::MatrixXi::RowXpr &Y_i,
                         boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_w = 0.0;
  Eigen::VectorXd w_s(2);
  Eigen::MatrixXd mlambda(2,G); // lambda for each cause
  Eigen::VectorXd cum_lambda(2); // cumulative lambda for each cause
  Eigen::MatrixXd tz(2,2); // row of z for each cause

  mlambda.row(0) = lambda0;
  mlambda.row(1) = lambda1;

  // std::cout << "Drawing...\n";
  for (int d = 0; d < 2; d++) {
      w_s(d) = stan::math::normal_rng(w(d), 1.0, rng);
  logr_w +=
      stan::math::normal_lpdf(w_s(d), 0.0, 1.0) -
      stan::math::normal_lpdf(w(d), 0.0, 1.0);
  }

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int k = 0; k < N; k++) {
    cum_lambda.setZero();
    tz.row(0) = z0.row(k);
    tz.row(1) = z1.row(k);

    for (int c = 0; c < 2; c++) {
      for (int g = 0; g < seg(k); g++) {
        cum_lambda(c) += len(g) * mlambda(c,g);
      }
      cum_lambda(c) += H(k) * mlambda(c,seg(k));

      logr_w -= cum_lambda(c) * ( stan::math::exp( beta(c) + theta(k,c) - gamma(c) * stan::math::distance(tz.row(c), w_s)) -
                                  stan::math::exp( beta(c) + theta(k,c) - gamma(c) * stan::math::distance(tz.row(c), w)));

      if (Y_i(k) == c) {
        logr_w -= gamma(c) * ( stan::math::distance(tz.row(c), w_s) - stan::math::distance(tz.row(c), w));
      }
    }
  }

  // accept or reject?
  if ((logr_w) > 0.0 || (logr_w >
                         stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
      acc_w += 1.0;
  }
  else {
      w_s = w;
  }
  return w_s;
}
