#ifndef __UPDATE_BETA_H_
#define __UPDATE_BETA_H_

#include <stan/math.hpp>

Eigen::VectorXd
update_beta(double &beta, double &acc_beta, const double &mu_beta,
            const double &sigma_beta, const double &jump_beta,
            const Eigen::MatrixXd::RowXpr &lambda,
            const Eigen::MatrixXd::ColXpr &theta, const double &gamma,
            const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
            const int &N, const Eigen::MatrixXi::RowXpr &NA,
            const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
            const Eigen::MatrixXd::RowXpr &H,
            const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
            boost::ecuyer1988 &rng);

// single beta update, return cumulative labmda!
Eigen::VectorXd update_single_beta(
    double &beta0, double &beta1, double &acc_beta, const double &mu_beta,
    const double &sigma_beta, const double &jump_beta,
    const Eigen::MatrixXd::RowXpr &lambda0,
    const Eigen::MatrixXd::RowXpr &lambda1, const Eigen::MatrixXd &theta,
    const Eigen::VectorXd &gamma, const Eigen::MatrixXd &z,
    const Eigen::MatrixXd::RowXpr &w0, const Eigen::MatrixXd::RowXpr &w1,
    const int &N, const Eigen::MatrixXi::RowXpr &NA, const Eigen::VectorXd &len,
    const Eigen::MatrixXi::RowXpr &seg, const Eigen::MatrixXd::RowXpr &H,
    boost::ecuyer1988 &rng);

Eigen::VectorXd update_cum_lambda(const Eigen::MatrixXd::RowXpr &lambda,
                                  const int &N,
                                  const Eigen::MatrixXi::RowXpr &NA,
                                  const Eigen::VectorXd &len,
                                  const Eigen::MatrixXi::RowXpr &seg,
                                  const Eigen::MatrixXd::RowXpr &H);

#endif // __UPDATE_BETA_H_
