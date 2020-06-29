#include <stan/math.hpp>
#include "update_gamma.h"

void update_gamma(Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
                  const Eigen::MatrixXd &lambda0, const Eigen::MatrixXd &lambda1,
                  const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta,
                  const Eigen::MatrixXd &z0, const Eigen::MatrixXd &z1,
                  const Eigen::MatrixXd &w,
                  const int &I, const int &N, const int &G, const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                  const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y,
                  boost::ecuyer1988 &rng) {
  
  // log-acceptance ratio
  Eigen::VectorXd logr_gamma(2); // gamma for each cause
  logr_gamma.setZero();
  Eigen::VectorXd gamma_s(2);
  Eigen::VectorXd cum_lambda(2); // cumulative lambda for each cause
  Eigen::MatrixXd lambda(I,G); 
  Eigen::MatrixXd z(N,2); 

  for (int c = 0; c < 2; c++) {
    if (c == 0) {
      lambda = lambda0;
      z = z0;
    }
    else {
      lambda = lambda1;
      z = z1;
    } 
    // std::cout << "Drawing...\n";
    gamma_s(c) = stan::math::lognormal_rng(stan::math::log(gamma(c)), 1.0, rng);

    // prior and jumping rule
    logr_gamma(c) +=
      stan::math::lognormal_lpdf(gamma_s(c), 0.0, 1.0) -
      stan::math::lognormal_lpdf(gamma(c), 0.0, 1.0) +
      stan::math::lognormal_lpdf(gamma_s(c), stan::math::log(gamma(c)), 1.0) -
      stan::math::lognormal_lpdf(gamma(c), stan::math::log(gamma_s(c)), 1.0);

    // std::cout << "Calculating the log-acceptance ratio...\n";
    for (int i = 0; i < I; i++) {
      for (int k = 0; k < N; k++) {
        cum_lambda.setZero();

        for (int g = 0; g < seg(i,k); g++) {
          cum_lambda(c) += len(g) * lambda(i,g);
        }
        cum_lambda(c) += H(i,k) * lambda(i,seg(i,k));

        logr_gamma(c) -= cum_lambda(c) * ( stan::math::exp( beta(i,c) + theta(k,c) - gamma_s(c) * stan::math::distance(z.row(k), w.row(i))) -
                                           stan::math::exp( beta(i,c) + theta(k,c) - gamma(c) * stan::math::distance(z.row(k), w.row(i))));

        if (Y(i,k) == c) {
          logr_gamma(c) += (gamma_s(c) - gamma(c)) * stan::math::distance(z.row(k), w.row(i));
        }
      }
    }

    // accept or reject?
    if ((logr_gamma(c)) > 0.0 || (logr_gamma(c) >
                                  stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
      acc_gamma(c) += 1.0;
      gamma(c) = gamma_s(c);
    }
  }
}
