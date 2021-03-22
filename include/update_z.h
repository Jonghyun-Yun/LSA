#ifndef __UPDATE_Z_H_
#define __UPDATE_Z_H_

#include <stan/math.hpp>

Eigen::VectorXd update_z(const Eigen::MatrixXd::RowXpr &z, double &acc_z,
                         const Eigen::MatrixXd::RowXpr &mu_z,
                         const Eigen::MatrixXd::RowXpr &sigma_z,
                         const Eigen::MatrixXd::RowXpr &jump_z,
                         const Eigen::VectorXd &cum_lambda,
                         const Eigen::MatrixXd::ColXpr &beta,
                         const double &theta, const double &gamma,
                         const Eigen::MatrixXd &w,
                         const int &I,
                         const Eigen::MatrixXi::ColXpr &NA,
                         const Eigen::VectorXd &len, const Eigen::MatrixXi::ColXpr &seg,
                         const Eigen::MatrixXd::ColXpr &H,
                         const Eigen::MatrixXi::ColXpr &Y_k, const int cause,
                         boost::ecuyer1988 &rng);

#endif // __UPDATE_Z_H_
