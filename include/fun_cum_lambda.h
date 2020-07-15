#ifndef __FUN_CUM_LAMBDA_H_
#define __FUN_CUM_LAMBDA_H_

void fun_cum_lambda(const Eigen::MatrixXd &lambda, Eigen::MatrixXd &cum_lambda,
                    const int &I, const int &N, const Eigen::MatrixXi &NA,
                    const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
                    const Eigen::MatrixXd &H);

#endif // __FUN_CUM_LAMBDA_H_
