#ifndef __UPDATE_Z_H_
#define __UPDATE_Z_H_

Eigen::VectorXd update_z(const Eigen::MatrixXd::RowXpr &z, double &acc_z,
                         const double &mu_z, const double &sigma_z, const double &jump_z,
                         const Eigen::MatrixXd &lambda,
                         const Eigen::MatrixXd::ColXpr &beta,
                         const double &theta, const double &gamma,
                         const Eigen::MatrixXd &w,
                         const int &I,
                         const Eigen::VectorXd &len, const Eigen::MatrixXi::ColXpr &seg,
                         const Eigen::MatrixXd::ColXpr &H,
                         const Eigen::MatrixXi::ColXpr &Y_k, const int cause,
                         boost::ecuyer1988 &rng);

#endif // __UPDATE_Z_H_
