#include <stan/math.hpp>
#include "par_fun_lp.h"

double fun_lp(const Eigen::MatrixXd &a_lambda, const Eigen::MatrixXd &b_lambda,
              const Eigen::MatrixXd &mu_beta, const Eigen::MatrixXd &sigma_beta,
              const Eigen::MatrixXd &mu_theta, const Eigen::MatrixXd &sigma_theta,
              const double &a_sigma, const double &b_sigma,
              const Eigen::VectorXd &mu_gamma, const Eigen::VectorXd &sigma_gamma,
              const Eigen::MatrixXd &mu_z, const Eigen::MatrixXd &sigma_z,
              const Eigen::MatrixXd &mu_w, const Eigen::MatrixXd &sigma_w,
              const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &cum_lambda,
              const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta, const double &sigma,
              const Eigen::VectorXd &gamma,
              const Eigen::MatrixXd &z, const Eigen::MatrixXd &w,
              const int &I, const int &N, const int &G,
              const Eigen::VectorXd &len, const Eigen::MatrixXi &seg,
              const Eigen::MatrixXd &H, const Eigen::MatrixXi &Y) {

    using stan::math::to_vector;
    using stan::math::distance;
    using stan::math::exp;
    using stan::math::square;
    using stan::math::normal_lpdf;
    using stan::math::lognormal_lpdf;
    using stan::math::gamma_lpdf;
    using stan::math::inv_gamma_lpdf;

    double lp_ = 0.0;
    // Eigen::MatrixXd cum_lambda(I,N); // cumulative lambda for each cause
    // Eigen::MatrixXd lambda(I,G);
    // Eigen::MatrixXd z(N,2);
    double log_prop_hazard; // to avoid double evaluation

    for (int c = 0; c < 2; c++) {

        // cum_lambda.setZero();
        // if (c == 0) {
        //     lambda = lambda0;
        //     z = z0;
        // }
        // else {
        //     lambda = lambda1;
        //     z = z1;
        // }
        
        // std::cout << "Calculating the log-acceptance ratio...\n";
        for (int i = 0; i < I; i++) {
            for (int k = 0; k < N; k++) {

                // for (int g = 0; g < seg(i,k); g++) {
                //     cum_lambda(i,k) += len(g) * lambda(i,g);
                // }
                // cum_lambda(i,k) += H(i,k) * lambda(i,seg(i,k));

                log_prop_hazard = beta(i,c) + theta(k,c) - gamma(c) * distance(z.row(c*N + k), w.row(i));
                lp_ -= cum_lambda(c*I + i,k) * exp( log_prop_hazard );

                if (Y(i,k) == c) {
                    lp_ += log(lambda(c * I + i, seg(i, k))) + log_prop_hazard;
                }
            }
        }
    }

        // priors...
    lp_ += inv_gamma_lpdf(square(sigma), a_sigma, b_sigma) +
        gamma_lpdf(to_vector(lambda.block(0,0,I,G)), to_vector(a_lambda), to_vector(b_lambda)) +
        gamma_lpdf(to_vector(lambda.block(I,0,I,G)), to_vector(a_lambda), to_vector(b_lambda)) +
        normal_lpdf(to_vector(beta), to_vector(mu_beta), to_vector(sigma_beta)) +
        normal_lpdf(to_vector(theta), to_vector(mu_theta), to_vector(sigma_theta)) +
        lognormal_lpdf(gamma, mu_gamma, sigma_gamma) +
        normal_lpdf(to_vector(z.block(0,0,N,2)), to_vector(mu_z), to_vector(sigma_z)) +
        normal_lpdf(to_vector(z.block(N,0,N,2)), to_vector(mu_z), to_vector(sigma_z)) +
        normal_lpdf(to_vector(w), to_vector(mu_w), to_vector(sigma_w));
 
    return lp_;
}



