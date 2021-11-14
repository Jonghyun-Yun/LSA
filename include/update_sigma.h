#ifndef __UPDATE_SIGMA_H_
#define __UPDATE_SIGMA_H_

#include <stan/math.hpp>

double update_sigma(const double& a_sigma, const double &b_sigma,
                           const Eigen::MatrixXd &theta, const int &N, boost::ecuyer1988 &rng);

#endif // __UPDATE_SIGMA_H_
