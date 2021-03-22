#include "update_beta.h"
#include <stan/math.hpp>

Eigen::VectorXd update_beta(double &beta, double &acc_beta, const double &mu_beta,
            const double &sigma_beta, const double &jump_beta,
            const Eigen::MatrixXd::RowXpr &lambda,
            const Eigen::MatrixXd::ColXpr &theta, const double &gamma,
            const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
            const int &N, const Eigen::MatrixXi::RowXpr &NA,
            const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
            const Eigen::MatrixXd::RowXpr &H,
            const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
            boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_beta = 0.0;
  double beta_s;
  Eigen::VectorXd cum_lambda = Eigen::VectorXd::Zero(N);

  // // std::cout << "Drawing...\n";
  beta_s = stan::math::normal_rng(beta, jump_beta, rng);
  // beta_s = 0;

  // prior and jump rule
  logr_beta += stan::math::normal_lpdf(beta_s, mu_beta, sigma_beta) -
               stan::math::normal_lpdf(beta, mu_beta, sigma_beta);

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
  for (int k = 0; k < N; k++) {

    if (NA(k) == 1) {

      for (int g = 0; g < seg(k); g++) {
        cum_lambda(k) += len(g) * lambda(g);
      }
      cum_lambda(k) += H(k) * lambda(seg(k));

      logr_beta -= cum_lambda(k) *
                   (stan::math::exp(beta_s + theta(k) -
                                    gamma * stan::math::distance(z.row(k), w)) -
                    stan::math::exp(beta + theta(k) -
                                    gamma * stan::math::distance(z.row(k), w)));

      if (Y_i(k) == cause) {
        logr_beta += beta_s - beta;
      }
    }
  }

  // accept or reject?
  if ((logr_beta > 0.0) ||
      (logr_beta > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    beta = beta_s;
    acc_beta += 1.0;
  }

  return cum_lambda;
}

// single beta update, return cumulative labmda!
Eigen::VectorXd update_single_beta(
    double &beta0, double &beta1, double &acc_beta,
    const double &mu_beta, const double &sigma_beta, const double &jump_beta,
    const Eigen::MatrixXd::RowXpr &lambda0,
    const Eigen::MatrixXd::RowXpr &lambda1, const Eigen::MatrixXd &theta,
    const Eigen::VectorXd &gamma, const Eigen::MatrixXd &z,
    const Eigen::MatrixXd::RowXpr &w0, const Eigen::MatrixXd::RowXpr &w1,
    const int &N, const Eigen::MatrixXi::RowXpr &NA, const Eigen::VectorXd &len,
    const Eigen::MatrixXi::RowXpr &seg, const Eigen::MatrixXd::RowXpr &H,
    boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_beta = 0.0;
  double beta_s;
  // Eigen::MatrixXd cum_lambda = Eigen::MatrixXd::Zero(2, N); // column-major storge!
  // https://eigen.tuxfamily.org/dox/group__TutorialReshapeSlicing.html
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> cum_lambda = Eigen::MatrixXd::Zero(2, N);
  Eigen::VectorXd res(2*N);

  Eigen::MatrixXd lambda(2, lambda0.cols()); // cumulative lambda for each cause
  Eigen::MatrixXd w(2, 2);      // w for each cause

  lambda.row(0) = lambda0;
  lambda.row(1) = lambda1;

  w.row(0) = w0;
  w.row(1) = w1;

  // // std::cout << "Drawing...\n";
  beta_s = stan::math::normal_rng(beta1, jump_beta, rng);
  // beta_s = 0.0;

  // prior and jump rule
  logr_beta += stan::math::normal_lpdf(beta_s, mu_beta, sigma_beta) -
               stan::math::normal_lpdf(beta1, mu_beta, sigma_beta); // prior

  // std::cout << "Calculating the log-acceptance ratio of lambda...\n";
    for (int k = 0; k < N; k++) {
      if (NA(k) == 1) {
        logr_beta += beta_s - beta1; // hazard

  for (int c = 0; c < 2; c++) {
        for (int g = 0; g < seg(k); g++) {
          cum_lambda(c, k) += len(g) * lambda(c, g);
        }
        cum_lambda(c, k) += H(k) * lambda(c, seg(k));

        logr_beta -=
            cum_lambda(c, k) *
            (stan::math::exp(
                 beta_s + theta(k, c) -
                 gamma(c) * stan::math::distance(z.row(c * N + k), w.row(c))) -
             stan::math::exp(
                 beta1 + theta(k, c) -
                 gamma(c) * stan::math::distance(z.row(c * N + k), w.row(c))));
      }
    }
  }

  // accept or reject?
  if ((logr_beta > 0.0) ||
      (logr_beta > stan::math::log(stan::math::uniform_rng(0.0, 1.0, rng)))) {
    beta0 = beta_s;
    beta1 = beta_s;
    acc_beta += 1.0;
  }
  res = Eigen::Map<Eigen::VectorXd>(cum_lambda.data(), 2 * N);
  // std::cout << lambda(0, 0) << std::endl
  //           << lambda(0, 1) << std::endl
  //           << cum_lambda(0, 0) << std::endl
  //           << res(0) << std::endl
  //           << "-----------------" << std::endl;
  return res;
  // return cum_lambda;
}

// udpate cum_lambda. no beta update
Eigen::VectorXd update_cum_lambda(
            const Eigen::MatrixXd::RowXpr &lambda,
            const int &N, const Eigen::MatrixXi::RowXpr &NA,
            const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
            const Eigen::MatrixXd::RowXpr &H) {

  Eigen::VectorXd cum_lambda = Eigen::VectorXd::Zero(N);

  for (int k = 0; k < N; k++) {

    if (NA(k) == 1) {

      for (int g = 0; g < seg(k); g++) {
        cum_lambda(k) += len(g) * lambda(g);
      }
      cum_lambda(k) += H(k) * lambda(seg(k));
    }
  }

  return cum_lambda;
}
