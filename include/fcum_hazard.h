#ifndef __FCUM_HAZARD_H_
#define __FCUM_HAZARD_H_

double fcum_hazard( const double &cum_lambda, const double &beta,
                    const double &theta,
                    const double &gamma, const Eigen::MatrixXd::RowXpr &z,
                    const Eigen::MatrixXd::RowXpr &w );

#endif // __FCUM_HAZARD_H_
