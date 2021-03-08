#include <stan/math.hpp>
#include "tbb/blocked_range2d.h"
#include "update_one_free_gamma.h"

void update_one_free_gamma(
    const int fc, Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
    const Eigen::VectorXd &mu_gamma, const Eigen::VectorXd &sigma_gamma,
    const Eigen::VectorXd &jump_gamma, const Eigen::MatrixXd &cum_lambda,
    const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta,
    const Eigen::MatrixXd &z, const Eigen::MatrixXd &w, const int &I,
    const int &N, const int &G, const Eigen::MatrixXi &NA,
    const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
    const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y, bool RUN_PAR,
    boost::ecuyer1988 &rng) {

  // fc: a "F"ree gamma variable response type "C", this is one to get updated.
  // do nothing for gamma(c) for c != fc
  int c = fc;

  // we draw a single free gamme gamma_s
  // Eigen::VectorXd gamma_s(2);
  double gamma_s;

  // log-acceptance ratio
  double logr_gamma = 0.0;

  // std::cout << "Drawing...\n";
  gamma_s =
      stan::math::lognormal_rng(stan::math::log(gamma(c)), jump_gamma(c), rng);

  // prior and jumping rule
  logr_gamma +=
      stan::math::lognormal_lpdf(gamma_s, mu_gamma(c), sigma_gamma(c)) -
      stan::math::lognormal_lpdf(gamma(c), mu_gamma(c), sigma_gamma(c)) +
      stan::math::lognormal_lpdf(gamma(c), stan::math::log(gamma_s),
                                 jump_gamma(c)) -
      stan::math::lognormal_lpdf(gamma_s, stan::math::log(gamma(c)),
                                 jump_gamma(c));

  if (RUN_PAR) {
    logr_gamma += tbb::parallel_reduce(
        tbb::blocked_range2d<int>(0, I, 0, N), 0.0,
        [&](tbb::blocked_range2d<int> r, double running_total) {
          for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
            for (int k = r.cols().begin(); k < r.cols().end(); ++k) {
              if (NA(i, k) == 1) {

                running_total -=
                    cum_lambda(c * I + i, k) *
                    (stan::math::exp(
                         beta(i, c) + theta(k, c) -
                         gamma_s * stan::math::distance(z.row(c * N + k),
                                                        w.row(c * I + i))) -
                     stan::math::exp(
                         beta(i, c) + theta(k, c) -
                         gamma(c) * stan::math::distance(z.row(c * N + k),
                                                         w.row(c * I + i))));

                if (Y(i, k) == c) {
                  running_total -=
                      (gamma_s - gamma(c)) *
                      stan::math::distance(z.row(c * N + k), w.row(c * I + i));
                }
              }
            }
          }
          return running_total;
        },
        std::plus<double>());

  } else {

    // for (int c = 0; c < 2; c++) {
    // std::cout << "Calculating the log-acceptance ratio...\n";
    for (int i = 0; i < I; i++) {
      for (int k = 0; k < N; k++) {

        if (NA(i, k) == 1) {

          // cum_lambda.setZero();

          // for (int g = 0; g < seg(i,k); g++) {
          //   cum_lambda(c) += len(g) * lambda(i,g);
          // }
          // cum_lambda(c) += H(i,k) * lambda(i,seg(i,k));

          logr_gamma -= cum_lambda(c * I + i, k) *
                        (stan::math::exp(
                             beta(i, c) + theta(k, c) -
                             gamma_s * stan::math::distance(z.row(c * N + k),
                                                            w.row(c * I + i))) -
                         stan::math::exp(beta(i, c) + theta(k, c) -
                                         gamma(c) * stan::math::distance(
                                                        z.row(c * N + k),
                                                        w.row(c * I + i))));

          if (Y(i, k) == c) {
            logr_gamma -=
                (gamma_s - gamma(c)) *
                stan::math::distance(z.row(c * N + k), w.row(c * I + i));
          }
        }
      }
    }
    // }
  }

  // accept or reject?
  if ((logr_gamma > 0.0) ||
      (logr_gamma > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    acc_gamma.array() += 1.0;
    gamma(c) = gamma_s;
  }
  // std::cout << "updating gamma...\n" << "logr_gamma: " << logr_gamma <<
  // "gamma and c" << gamma(c) << ", " << c << std::endl;
}
