#include <stan/math.hpp>
#include "tbb/blocked_range3d.h"
#include "par_fun_lp.h"

double par_fun_lp(const Eigen::MatrixXd &a_lambda, const Eigen::MatrixXd &b_lambda,
                  const Eigen::MatrixXd &mu_beta, const Eigen::MatrixXd &sigma_beta,
                  const Eigen::MatrixXd &mu_theta, const Eigen::MatrixXd &sigma_theta,
                  const double &a_sigma, const double &b_sigma,
                  const Eigen::VectorXd &mu_gamma, const Eigen::VectorXd &sigma_gamma,
                  const Eigen::MatrixXd &mu_z, const Eigen::MatrixXd &sigma_z,
                  const Eigen::MatrixXd &mu_w, const Eigen::MatrixXd &sigma_w,
                  const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &cum_lambda,
                  const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta, const double &sigma,
                  const Eigen::VectorXd &gamma,
                  const Eigen::MatrixXd &z, const Eigen::MatrixXd &w,
                  const int &I, const int &N, const int &G,
                  const Eigen::MatrixXi &NA,
                  const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                  const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y,
                  bool SINGLE_Z, bool SINGLE_W, bool UPDATE_GAMMA, int ONE_FREE_GAMMA) {

    using stan::math::to_vector;
    using stan::math::distance;
    using stan::math::exp;
    using stan::math::log;
    using stan::math::square;
    using stan::math::normal_lpdf;
    using stan::math::lognormal_lpdf;
    using stan::math::gamma_lpdf;
    using stan::math::inv_gamma_lpdf;

    double lp_ = 0.0;
    // Eigen::MatrixXd cum_lambda(I,N); // cumulative lambda for each cause
    // Eigen::MatrixXd lambda(I,G);
    // Eigen::MatrixXd z(N,2);
    // double log_prop_hazard; // to avoid double evaluation

    lp_ += tbb::parallel_reduce(
        tbb::blocked_range3d<int>(0, 2, 0, I, 0, N), 0.0,
        [&](tbb::blocked_range3d<int> r, double running_total) {
          for (int c = r.pages().begin(); c < r.pages().end(); ++c) {
            for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
              for (int k = r.cols().begin(); k < r.cols().end(); ++k) {
                // for (int g = 0; g < seg(i,k); g++) {
                //     cum_lambda(c*I + i,k) += len(g) * lambda(i,g);
                // }
                // cum_lambda(i,k) += H(i,k) * lambda(i,seg(i,k));

                if (NA(i, k) == 1) {

                  running_total -= cum_lambda(c * I + i, k) *
                                   exp(beta(i, c) + theta(k, c) -
                                       gamma(c) * distance(z.row(c * N + k),
                                                           w.row(c * I + i)));

                  if (Y(i, k) == c) {
                    running_total +=
                        log(lambda(c * I + i, seg(i, k))) + beta(i, c) +
                        theta(k, c) -
                        gamma(c) * distance(z.row(c * N + k), w.row(c * I + i));
                  }
                }
              }
            }
          }
          return running_total;
        },
        std::plus<double>());

    // priors...
    lp_ += inv_gamma_lpdf(square(sigma), a_sigma, b_sigma) +
        gamma_lpdf(to_vector(lambda.block(0,0,I,G)), to_vector(a_lambda), to_vector(b_lambda)) +
        gamma_lpdf(to_vector(lambda.block(I,0,I,G)), to_vector(a_lambda), to_vector(b_lambda)) +
        normal_lpdf(to_vector(beta), to_vector(mu_beta), to_vector(sigma_beta)) +
        normal_lpdf(to_vector(theta), to_vector(mu_theta), to_vector(sigma_theta));

    if (SINGLE_Z) {
        lp_ += normal_lpdf(to_vector(z.block(0,0,N,2)), to_vector(mu_z), to_vector(sigma_z));
    }
    else {
        lp_ += normal_lpdf(to_vector(z.block(0,0,N,2)), to_vector(mu_z), to_vector(sigma_z)) +
        normal_lpdf(to_vector(z.block(N,0,N,2)), to_vector(mu_z), to_vector(sigma_z));
    }

    if (SINGLE_W) {
       lp_ += normal_lpdf(to_vector(w.block(0,0,I,2)), to_vector(mu_w), to_vector(sigma_w));
    }
    else {
        lp_ += normal_lpdf(to_vector(w.block(0,0,I,2)), to_vector(mu_w), to_vector(sigma_w)) +
            normal_lpdf(to_vector(w.block(I,0,I,2)), to_vector(mu_w), to_vector(sigma_w));
    }

    if (UPDATE_GAMMA) {
        if (SINGLE_Z && SINGLE_W) {
            // abs to play with sign
            if ((ONE_FREE_GAMMA == 99) || (ONE_FREE_GAMMA == 51)) {
                lp_ += lognormal_lpdf(std::abs(gamma(1)), mu_gamma(1), sigma_gamma(1));
            }
            else {
                lp_ += lognormal_lpdf(std::abs(gamma(ONE_FREE_GAMMA)),
                                      mu_gamma(ONE_FREE_GAMMA), sigma_gamma(ONE_FREE_GAMMA));
            }
        } else if (ONE_FREE_GAMMA == 51) {
            lp_ += lognormal_lpdf(std::abs(gamma(1)), mu_gamma(1), sigma_gamma(1));
        }
        else {
            lp_ += lognormal_lpdf(gamma.array().abs(), mu_gamma, sigma_gamma); // gamma can be negative.
        }
    }

    return lp_;
}



