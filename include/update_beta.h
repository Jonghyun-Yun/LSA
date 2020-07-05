#ifndef __UPDATE_BETA_H_
#define __UPDATE_BETA_H_

Eigen::VectorXd update_beta(double &beta, double &acc_beta,
                            const double &mu_beta, const double &sigma_beta, const double &jump_beta,
                            const Eigen::MatrixXd::RowXpr &lambda,
                            const Eigen::MatrixXd::ColXpr &theta,
                            const double &gamma,
                            const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
                            const int &N,
                            const Eigen::VectorXd &len, const Eigen::MatrixXi::RowXpr &seg,
                            const Eigen::MatrixXd::RowXpr &H,
                            const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                            boost::ecuyer1988 &rng);

#endif // __UPDATE_BETA_H_
