#include <stan/math.hpp>
#include "fcum_hazard.h"

// this function cannot be used if z or w is passed from another function calls this one
// https://stackoverflow.com/questions/48856518/eigen-library-compiler-error-when-passing-block-reference-to-templated-function
// if you want to use this one, you need to bind z or w to a vector.
double fcum_hazard( const double &cum_lambda, const double &beta,
                    const double &theta,
                    const double &gamma, const Eigen::MatrixXd::RowXpr &z,
                    const Eigen::MatrixXd::RowXpr &w ) {

    return (cum_lambda * stan::math::exp( beta + theta - gamma * stan::math::distance(z, w) ));

}
