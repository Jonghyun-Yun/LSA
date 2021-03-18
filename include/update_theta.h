#ifndef __UPDATE_THETA_H_
#define __UPDATE_THETA_H_

#include <stan/math.hpp>

void update_theta(double &theta, double &acc_theta,
                  const double &mu_theta, const double &jump_theta,
                  const double &sigma,
                  const Eigen::VectorXd &cum_lambda,
                  const Eigen::MatrixXd::ColXpr &beta,
                  const double &gamma,
                  const Eigen::MatrixXd::RowXpr &z, const Eigen::MatrixXd &w,
                  const int &I,
                  const Eigen::MatrixXi::ColXpr &NA,
                  const Eigen::VectorXd &len, const Eigen::MatrixXi::ColXpr &seg,
                  const Eigen::MatrixXd::ColXpr &H,
                  const Eigen::MatrixXi::ColXpr &Y_k, const int cause,
                  boost::ecuyer1988 &rng);

#endif // __UPDATE_THETA_H_
