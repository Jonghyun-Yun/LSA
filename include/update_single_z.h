#ifndef __UPDATE_SINGLE_Z_H_
#define __UPDATE_SINGLE_Z_H_

Eigen::VectorXd update_single_z(const Eigen::MatrixXd::RowXpr &z, double &acc_z,
                                const double &mu_z, const double &sigma_z, const double &jump_z,
                                const Eigen::MatrixXd::ColXpr &cum_lambda,
                                const Eigen::MatrixXd &beta,
                                const Eigen::MatrixXd::RowXpr &theta,
                                const Eigen::VectorXd &gamma,
                                const Eigen::MatrixXd &w,
                                const int &I,
                                const Eigen::MatrixXi::ColXpr &NA,
                                const Eigen::VectorXd &len, const Eigen::MatrixXi::ColXpr &seg,
                                const Eigen::MatrixXd::ColXpr &H,
                                const Eigen::MatrixXi::ColXpr &Y_k,
                                boost::ecuyer1988 &rng);

#endif // __UPDATE_SINGLE_Z_H_
