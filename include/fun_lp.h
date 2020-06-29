#ifndef __FUN_LP_H_
#define __FUN_LP_H_

double fun_lp(const Eigen::MatrixXd &a_lambda, const Eigen::MatrixXd &b_lambda,
              const Eigen::MatrixXd &mu_beta, const Eigen::MatrixXd &sigma_beta,
              const Eigen::MatrixXd &mu_theta, const Eigen::MatrixXd &sigma_theta,
              const double &a_sigma, const double &b_sigma,
              const Eigen::VectorXd &mu_gamma, const Eigen::VectorXd &sigma_gamma,
              const Eigen::MatrixXd &mu_z, const Eigen::MatrixXd &sigma_z,
              const Eigen::MatrixXd &mu_w, const Eigen::MatrixXd &sigma_w,
              const Eigen::MatrixXd &lambda0, const Eigen::MatrixXd &lambda1,
              const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta, const double &sigma,
              const Eigen::VectorXd &gamma,
              const Eigen::MatrixXd &z0, const Eigen::MatrixXd &z1, const Eigen::MatrixXd &w,
              const int &I, const int &N, const int &G,
              const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
              const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y);

#endif // __FUN_LP_H_
