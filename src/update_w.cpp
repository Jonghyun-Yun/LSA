#include <stan/math.hpp>
#include "update_w.h"

Eigen::VectorXd update_w(const Eigen::MatrixXd::RowXpr &w, double &acc_w,
                         const Eigen::MatrixXd::RowXpr &lambda,
                         const Eigen::MatrixXd::ColXpr &theta,
                         const double &gamma,
                         const Eigen::MatrixXd &z,
                         const int &N,
                         const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                         const Eigen::MatrixXd::RowXpr &H,
                         const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                         boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_w = 0.0;
  Eigen::VectorXd w_s(2);
  double cum_lambda;

  // std::cout << "Drawing...\n";
  for (int d = 0; d < 2; d++) {
      z_s(d) = stan::math::normal_rng(z(d), 1.0, rng);
  logr_z +=
      stan::math::normal_lpdf(z_s(d), 0.0, 1.0) -
      stan::math::normal_lpdf(z(d), 0.0, 1.0);
  }

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int i = 0; i < I; i++) {
  cum_lambda = 0.0;
    for (int g = 0; g < seg(i); g++) {
      cum_lambda += len(g) * lambda(i,g);
    }
    cum_lambda += H(i) * lambda(i,seg(i));

    logr_z -= cum_lambda * ( stan::math::exp( beta(i) + theta - gamma * stan::math::distance(z, w.row(i))) -
                             stan::math::exp( beta(i) + theta - gamma * stan::math::distance(z_s, w.row(i))));

      if (Y_k(i) == cause) {
          logr_z -= gamma * ( stan::math::distance(z, w.row(i)) - stan::math::distance(z_s, w.row(i)));
      }
    }

  // accept or reject?
  if ((logr_z) > 0.0 || (logr_z >
                              stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
      acc_z += 1.0;
      return z_s;
  }
  else {
      z_s = z;
      return z_s;
  }
}
