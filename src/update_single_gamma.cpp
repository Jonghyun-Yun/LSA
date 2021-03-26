#include <stan/math.hpp>
#include "tbb/blocked_range3d.h"
#include "update_single_gamma.h"

void update_single_gamma(Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
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
                         bool RUN_PAR,
                         boost::ecuyer1988 &rng) {


// for single gamma, gamma(1) > 0 and gamma(0) = -1.0 * gamma(1);
// we draw gamma_s(1), and set gamma_s(0) = -1.0 * gamma_s(1);
    Eigen::VectorXd gamma_s(2);

    // log-acceptance ratio
    double logr_gamma = 0.0;

    double gamma_1; // when gamma(1) < 0

    if (gamma(1) > 0) {
    // std::cout << "Drawing...\n";
    gamma_s(1) = stan::math::lognormal_rng(stan::math::log(gamma(1)), jump_gamma(1), rng);

    // prior and jumping rule
    logr_gamma +=
        stan::math::lognormal_lpdf(gamma_s(1), mu_gamma(1), sigma_gamma(1)) -
        stan::math::lognormal_lpdf(gamma(1), mu_gamma(1), sigma_gamma(1)) +
        stan::math::lognormal_lpdf(gamma(1), stan::math::log(gamma_s(1)), jump_gamma(1)) -
        stan::math::lognormal_lpdf(gamma_s(1), stan::math::log(gamma(1)), jump_gamma(1));
    } else {
       gamma_1 = std::abs(gamma(1));
       gamma_s(1) = stan::math::lognormal_rng(stan::math::log(gamma_1), jump_gamma(1), rng);
    // prior and jumping rule
    logr_gamma +=
        stan::math::lognormal_lpdf(gamma_s(1), mu_gamma(1), sigma_gamma(1)) -
        stan::math::lognormal_lpdf(gamma_1, mu_gamma(1), sigma_gamma(1)) +
        stan::math::lognormal_lpdf(gamma_1, stan::math::log(gamma_s(1)), jump_gamma(1)) -
        stan::math::lognormal_lpdf(gamma_s(1), stan::math::log(gamma_1), jump_gamma(1));

        gamma_s(1) *= -1.0; // put it back to negative
    }

    // respect gamma(0) sign
    if (gamma(0) < 0) {
      gamma_s(0) = -1.0 * std::abs(gamma_s(1));
    } else {
      gamma_s(0) = std::abs(gamma_s(1));
    }

    if (RUN_PAR) {
        
        logr_gamma += tbb::parallel_reduce(
                tbb::blocked_range3d<int>(0, 2, 0, I, 0, N),
                0.0,
                [&](tbb::blocked_range3d<int> r, double running_total)
                {
                    for(int c=r.pages().begin(); c < r.pages().end(); ++c)
                    {
                        for(int i=r.rows().begin(); i<r.rows().end(); ++i)
                        {
                            for(int k=r.cols().begin(); k<r.cols().end(); ++k)
                            {
                                if (NA(i,k) == 1) {

                                    running_total -=
                                        cum_lambda(c * I + i, k) *
                                        (stan::math::exp(beta(i, c) + theta(k, c) -
                                                         gamma_s(c) *
                                                         stan::math::distance(z.row(c*N + k), w.row(c*I + i))) -
                                         stan::math::exp(beta(i, c) + theta(k, c) -
                                                         gamma(c) *
                                                         stan::math::distance(z.row(c*N + k), w.row(c*I + i))));

                                    if (Y(i, k) == c) {
                                        running_total -= (gamma_s(c) - gamma(c)) *
                                            stan::math::distance(z.row(c*N + k), w.row(c*I + i));
                                    }
                                }
                            }
                        }
                    }
                    return running_total;
                }, std::plus<double>() );

    }
    else {

    for (int c = 0; c < 2; c++) {
        // std::cout << "Calculating the log-acceptance ratio...\n";
        for (int i = 0; i < I; i++) {
            for (int k = 0; k < N; k++) {

                if (NA(i,k) == 1) {

                    // cum_lambda.setZero();

                    // for (int g = 0; g < seg(i,k); g++) {
                    //   cum_lambda(c) += len(g) * lambda(i,g);
                    // }
                    // cum_lambda(c) += H(i,k) * lambda(i,seg(i,k));

                    logr_gamma -=
                        cum_lambda(c * I + i, k) *
                        (stan::math::exp(beta(i, c) + theta(k, c) -
                                         gamma_s(c) *
                                         stan::math::distance(z.row(c*N + k), w.row(c*I + i))) -
                         stan::math::exp(beta(i, c) + theta(k, c) -
                                         gamma(c) *
                                         stan::math::distance(z.row(c*N + k), w.row(c*I + i))));

                    if (Y(i, k) == c) {
                        logr_gamma -= (gamma_s(c) - gamma(c)) * stan::math::distance(z.row(c*N + k), w.row(c*I + i));
                    }
                }
            }
        }
    }
    }

    // accept or reject?
    if ((logr_gamma > 0.0) || (logr_gamma >
                               stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
        acc_gamma.array() += 1.0;
        gamma = gamma_s;
    }
}
