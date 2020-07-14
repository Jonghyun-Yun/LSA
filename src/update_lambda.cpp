#include <stan/math.hpp>
#include "update_lambda.h"

double update_lambda(const double &a_lambda,
                     const double &b_lambda,
                     const int g, const double &beta,
                     const Eigen::MatrixXd::ColXpr &theta, const double &gamma,
                     const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
                     const int &N, const Eigen::MatrixXi::RowXpr &NA,
                     const double &len, const Eigen::MatrixXi::RowXpr &seg,
                     const Eigen::MatrixXd::RowXpr &H,
                     const double &IY,
                     boost::ecuyer1988 &rng) {

    double post_a = a_lambda + 1.0 + IY;
    double post_b = b_lambda;

    for (int k = 0; k < N; k++) {

      if (NA(k) == 1) {

        if (g < seg(k)) {
          post_b +=
              len * std::exp(beta + theta(k) -
                              gamma * stan::math::distance(z.row(k), w));
        } else if (g == seg(k)) {
          post_b +=
              H(k) * std::exp(beta + theta(k) -
                              gamma * stan::math::distance(z.row(k), w));
        }
      }
    }

   return stan::math::gamma_rng(post_a, post_b, rng);

}
