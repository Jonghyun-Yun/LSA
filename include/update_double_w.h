#ifndef __UPDATE_DOUBLE_W_H_
#define __UPDATE_DOUBLE_W_H_

#include <stan/math.hpp>

Eigen::VectorXd update_double_w(const Eigen::MatrixXd::RowXpr &w, double &acc_w,
                                const Eigen::MatrixXd::RowXpr &mu_w,
                                const Eigen::MatrixXd::RowXpr &sigma_w,
                                const Eigen::MatrixXd::RowXpr &jump_w,
                                const Eigen::VectorXd &cum_lambda,
                                const double &beta,
                                const Eigen::MatrixXd::ColXpr &theta,
                                const double &gamma,
                                const Eigen::MatrixXd &z,
                                const int &N,
                                const Eigen::MatrixXi::RowXpr &NA,
                                const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                                const Eigen::MatrixXd::RowXpr &H,
                                const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                                boost::ecuyer1988 &rng);

#endif // __UPDATE_DOUBLE_W_H_
