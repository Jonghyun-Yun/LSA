#include <stan/math.hpp>
#include "update_single_z.h"

Eigen::VectorXd update_single_z(const Eigen::MatrixXd::RowXpr &z, double &acc_z,
                                const double &mu_z, const double &sigma_z, const double &jump_z,
                                const Eigen::MatrixXd::ColXpr &cum_lambda,
                                const Eigen::MatrixXd &beta,
                                const Eigen::MatrixXd::RowXpr &theta,
                                const Eigen::VectorXd &gamma,
                                const Eigen::MatrixXd &w,
                                const int &I,
                                const Eigen::MatrixXi::ColXpr &NA,
                                const Eigen::VectorXd &len, const Eigen::MatrixXi::ColXpr &seg,
                                const Eigen::MatrixXd::ColXpr &H,
                                const Eigen::MatrixXi::ColXpr &Y_k,
                                boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_z = 0.0;
  Eigen::VectorXd z_s(2);

  // std::cout << "Drawing...\n";
  for (int d = 0; d < 2; d++) {
    z_s(d) = stan::math::normal_rng(z(d), jump_z, rng);
    logr_z +=
      stan::math::normal_lpdf(z_s(d), mu_z, sigma_z) -
      stan::math::normal_lpdf(z(d), mu_z, sigma_z);
  }

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int i = 0; i < I; i++) {

    if(NA(i) == 1) {

      for (int c = 0; c < 2; c++) {
        // cum_lambda = 0.0;
        //   for (int g = 0; g < seg(i); g++) {
        //     cum_lambda += len(g) * lambda(i,g);
        //   }
        //   cum_lambda += H(i) * lambda(i,seg(i));

        logr_z -=
            cum_lambda(c * I + i) *
            (stan::math::exp(beta(i, c) + theta(c) -
                             gamma(c) * stan::math::distance(z_s, w.row(i))) -
             stan::math::exp(beta(i, c) + theta(c) -
                             gamma(c) * stan::math::distance(z, w.row(i))));

        if (Y_k(i) == c) {
          logr_z -= gamma(c) * ( stan::math::distance(z_s, w.row(i)) - stan::math::distance(z, w.row(i)));
        }
      }
    }
  }

  // accept or reject?
  if ((logr_z > 0.0) ||
      (logr_z > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    acc_z += 1.0;
  }
  else {
    z_s = z;
  }
  return z_s; // cannot assign values to z
}
