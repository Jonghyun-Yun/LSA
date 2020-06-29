#ifndef __UPDATE_GAMMA_H_
#define __UPDATE_GAMMA_H_

void update_gamma(Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
                  const Eigen::MatrixXd &lambda0, const Eigen::MatrixXd &lambda1,
                  const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta,
                  const Eigen::MatrixXd &z0, const Eigen::MatrixXd &z1, const Eigen::MatrixXd &w,
                  const int &I, const int &N, const int &G,
                  const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                  const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y,
                  boost::ecuyer1988 &rng);

#endif // __UPDATE_GAMMA_H_
