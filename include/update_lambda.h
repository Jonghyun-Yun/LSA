#ifndef __UPDATE_LAMBDA_H_
#define __UPDATE_LAMBDA_H_

void update_lambda(double &lambda, double &acc_lambda, const int g,
                   const double &beta, const Eigen::MatrixXd::ColXpr &theta,
                   const double &gamma,
                   const Eigen::MatrixXd &z, const Eigen::MatrixXd::RowXpr &w,
                   const int &N,
                   const double &len, const Eigen::MatrixXi::RowXpr &seg,
                   const Eigen::MatrixXd::RowXpr &H,
                   const Eigen::MatrixXi::RowXpr &Y_i, const int cause,
                   boost::ecuyer1988 &rng);

#endif // __UPDATE_LAMBDA_H_
