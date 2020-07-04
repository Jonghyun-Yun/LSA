// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;
// using Numer::integrate;
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using std::log;
using std::exp;

// typedef Eigen::Map<Eigen::MatrixXd> MapMat;
// typedef Eigen::Map<Eigen::VectorXd> MapVec;

// P(0.3 < X < 0.8), X ~ Beta(a, b)
class incirate: public Func
{
    private:
        bool T_F;
        MatrixXd lambda;
        VectorXd beta;
        VectorXd theta;
        VectorXd gamma;
        MatrixXd z;
        VectorXd w;
        VectorXd sj;
        VectorXd H;
        int G;
    public:
        incirate(bool T_F_, MatrixXd lambda_,
                 VectorXd beta_, VectorXd theta_, VectorXd gamma_,
                 MatrixXd z_, VectorXd w_, VectorXd sj_, VectorXd H_) :
            T_F(T_F_), lambda(lambda_), beta(beta_), theta(theta_),
            gamma(gamma_), z(z_), w(w_), sj(sj_), H(H_){
            G = H.size();
        }

        double operator()(const double& t) const {
            VectorXd cum_lambda(2);
            cum_lambda.setZero();

            double res = 1.0;

            int g = 0;
            while (t > sj(g+1)) {
                cum_lambda += H(g) * lambda.col(g);
                g++;
            }
            for (int c=0; c<2; c++) {
                res *= exp( -1.0 * ( cum_lambda(c) + ( t - sj(g) ) * lambda(c,g) ) * exp( beta(c) + theta(c) - gamma(c) * (z.row(c)-w).norm() ));
            }
            if (T_F) {
                res *= lambda(1,g) * exp( beta(1) + theta(1) - gamma(1) * (z.row(1)-w).norm() );
            }
            else {
                res *= lambda(0,g) * exp( beta(0) + theta(0) - gamma(0) * (z.row(0)-w).norm() );
            }
            return res;
        }

        double hazard(const double &t, const bool T_F) {
            double res;
            int g = 0;
            while (t > sj(g)) {
                g++;
            }
            if (T_F) {
                res = lambda(1,g) * exp( beta(1) + theta(1) - gamma(1) * (z.row(1)-w).norm() );
            }
            else {
                res = lambda(0,g) * exp( beta(0) + theta(0) - gamma(0) * (z.row(0)-w).norm() );
            }
            return res;
        }

        double accuracy(const double &t) {
            return hazard(t, 1) / ( hazard(t, 1) + hazard(t, 0) );
        }
};

// [[Rcpp::export]]
Rcpp::List integrate_incirate(double lower, double upper, bool T_F, MatrixXd lambda,
                              VectorXd beta,  VectorXd theta, VectorXd gamma,
                              MatrixXd z, VectorXd w, VectorXd sj, VectorXd H)
{
    incirate f(T_F, lambda, beta, theta, gamma, z, w, sj, H);
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return Rcpp::List::create(
        Rcpp::Named("approximate") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}

// [[Rcpp::export]]
VectorXd eval_incirate(VectorXd t, bool T_F, MatrixXd lambda,
                    VectorXd beta,  VectorXd theta, VectorXd gamma,
                    MatrixXd z, VectorXd w, VectorXd sj, VectorXd H)
{
    incirate f(T_F, lambda, beta, theta, gamma, z, w, sj, H);
    VectorXd res(t.size());
    for (int i=0; i<t.size(); i++) {
        res(i) = f(t(i));
    }
    return res;
}

// [[Rcpp::export]]
VectorXd eval_accuracy(VectorXd t, MatrixXd lambda,
                    VectorXd beta,  VectorXd theta, VectorXd gamma,
                    MatrixXd z, VectorXd w, VectorXd sj, VectorXd H)
{
    incirate f(1, lambda, beta, theta, gamma, z, w, sj, H);
    VectorXd res(t.size());
    for (int i=0; i<t.size(); i++) {
        res(i) = f.accuracy(t(i));
    }
    return res;
}
