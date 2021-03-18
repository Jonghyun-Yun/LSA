#include "update_w.h"

Eigen::VectorXd
update_w(const Eigen::MatrixXd::RowXpr &w, double &acc_w,
         const Eigen::MatrixXd::RowXpr &mu_w,
         const Eigen::MatrixXd::RowXpr &sigma_w,
         const Eigen::MatrixXd::RowXpr &jump_w,
         const Eigen::MatrixXd::RowXpr &cum_lambda0,
         const Eigen::MatrixXd::RowXpr &cum_lambda1,
         const Eigen::MatrixXd::RowXpr &beta, const Eigen::MatrixXd &theta,
         const Eigen::VectorXd &gamma, const Eigen::MatrixXd &z, const int &N,
         const Eigen::MatrixXi::RowXpr &NA, const Eigen::VectorXd &len,
         const Eigen::MatrixXi::RowXpr &seg, const Eigen::MatrixXd::RowXpr &H,
         const Eigen::MatrixXi::RowXpr &Y_i, boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_w = 0.0;
  Eigen::VectorXd w_s(2);
  // Eigen::MatrixXd lambda(2,G); // lambda for each cause
  Eigen::MatrixXd cum_lambda(2, N); // cumulative lambda for each cause
  // Eigen::MatrixXd tz(2,2); // row of z for each cause

  cum_lambda.row(0) = cum_lambda0;
  cum_lambda.row(1) = cum_lambda1;

  // std::cout << "Drawing...\n";
  for (int d = 0; d < 2; d++) {
    // symmetric jumping
    w_s(d) = stan::math::normal_rng(w(d), jump_w(d), rng);
    // prior
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
      // cum_lambda.setZero();
      // tz.row(0) = z0.row(k);
      // tz.row(1) = z1.row(k);

      // for single_z z.row(k) = z.row(N + k)
      for (int c = 0; c < 2; c++) {
        // for (int g = 0; g < seg(k); g++) {
        // cum_lambda(c) += len(g) * lambda(c,g);
        // }
        // cum_lambda(c) += H(k) * lambda(c,seg(k));

        logr_w -= cum_lambda(c, k) *
                  (stan::math::exp(
                       beta(c) + theta(k, c) -
                       gamma(c) * stan::math::distance(z.row(c * N + k), w_s)) -
                   stan::math::exp(
                       beta(c) + theta(k, c) -
                       gamma(c) * stan::math::distance(z.row(c * N + k), w)));

        if (Y_i(k) == c) {
          logr_w -= gamma(c) * (stan::math::distance(z.row(c * N + k), w_s) -
                                stan::math::distance(z.row(c * N + k), w));
        }
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

// update w_bgein (draw from a unit circle centered at z_begin)
// single z is assumed
Eigen::VectorXd update_w_begin(
    const Eigen::MatrixXd::RowXpr &w, double &acc_w,
    const Eigen::MatrixXd::RowXpr &mu_w, const Eigen::MatrixXd::RowXpr &sigma_w,
    const Eigen::MatrixXd::RowXpr &jump_w,
    const Eigen::MatrixXd::RowXpr &cum_lambda0,
    const Eigen::MatrixXd::RowXpr &cum_lambda1,
    const Eigen::MatrixXd::RowXpr &beta, const Eigen::MatrixXd &theta,
    const Eigen::VectorXd &gamma, const Eigen::MatrixXd &z, const int &N,
    const Eigen::MatrixXi::RowXpr &NA, const Eigen::VectorXd &len,
    const Eigen::MatrixXi::RowXpr &seg, const Eigen::MatrixXd::RowXpr &H,
    const Eigen::MatrixXi::RowXpr &Y_i, boost::ecuyer1988 &rng) {

  // log-acceptance ratio
  double logr_w = 0.0;
  Eigen::VectorXd w_s(2);
  // Eigen::MatrixXd lambda(2,G); // lambda for each cause
  Eigen::MatrixXd cum_lambda(2, N); // cumulative lambda for each cause
  // Eigen::MatrixXd tz(2,2); // row of z for each cause

  cum_lambda.row(0) = cum_lambda0;
  cum_lambda.row(1) = cum_lambda1;

  // std::cout << "Drawing...\n";
  // degree to sample
  double rad = stan::math::uniform_rng(0.0, 2.0 * stan::math::pi(), rng);
  w_s(0) = stan::math::cos(rad) + z(0, 0);
  w_s(1) = stan::math::sin(rad) + z(0, 1);


  // std::cout << "============================" << std::endl;
  // std::cout << "w: " << w_s(0) << "," << w_s(1) << std::endl;
  // std::cout << "z: " << z(0,0) << "," << z(0,1) << std::endl;

  for (int d = 0; d < 2; d++) {
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
      // cum_lambda.setZero();
      // tz.row(0) = z0.row(k);
      // tz.row(1) = z1.row(k);

      // for single_z z.row(k) = z.row(N + k)
      for (int c = 0; c < 2; c++) {
        // for (int g = 0; g < seg(k); g++) {
        // cum_lambda(c) += len(g) * lambda(c,g);
        // }
        // cum_lambda(c) += H(k) * lambda(c,seg(k));

        logr_w -= cum_lambda(c, k) *
                  (stan::math::exp(
                       beta(c) + theta(k, c) -
                       gamma(c) * stan::math::distance(z.row(c * N + k), w_s)) -
                   stan::math::exp(
                       beta(c) + theta(k, c) -
                       gamma(c) * stan::math::distance(z.row(c * N + k), w)));

        if (Y_i(k) == c) {
          logr_w -= gamma(c) * (stan::math::distance(z.row(c * N + k), w_s) -
                                stan::math::distance(z.row(c * N + k), w));
        }
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
