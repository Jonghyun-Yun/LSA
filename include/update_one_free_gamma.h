#ifndef __UPDATE_ONE_FREE_GAMMA_H_
#define __UPDATE_ONE_FREE_GAMMA_H_

void update_one_free_gamma(const int fc, Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
                         const Eigen::VectorXd &mu_gamma,
                         const Eigen::VectorXd &sigma_gamma,
                         const Eigen::VectorXd &jump_gamma,
                         const Eigen::MatrixXd &cum_lambda,
                         const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta,
                         const Eigen::MatrixXd &z, const Eigen::MatrixXd &w,
                         const int &I, const int &N, const int &G,
                         const Eigen::MatrixXi &NA,
                         const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                         const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y,
                         bool RUN_PAR,
                         boost::ecuyer1988 &rng);


#endif
