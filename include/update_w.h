#ifndef __UPDATE_W_H_
#define __UPDATE_W_H_

#include <stan/math.hpp>

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
         const Eigen::MatrixXi::RowXpr &Y_i, boost::ecuyer1988 &rng);

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
    const Eigen::MatrixXi::RowXpr &Y_i, boost::ecuyer1988 &rng);

#endif // __UPDATE_W_H_
