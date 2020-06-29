#ifndef __UPDATE_W_H_
#define __UPDATE_W_H_

Eigen::VectorXd update_w(const Eigen::MatrixXd::RowXpr &w, double &acc_w,
                         const Eigen::MatrixXd::RowXpr &lambda0,
                         const Eigen::MatrixXd::RowXpr &lambda1,
                         const Eigen::MatrixXd::RowXpr &beta,
                         const Eigen::MatrixXd &theta,
                         const Eigen::VectorXd &gamma,
                         const Eigen::MatrixXd &z0, const Eigen::MatrixXd &z1,
                         const int &N, const int &G, const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                         const Eigen::MatrixXd::RowXpr &H, const Eigen::MatrixXi::RowXpr &Y_i,
                         boost::ecuyer1988 &rng);

#endif // __UPDATE_W_H_
