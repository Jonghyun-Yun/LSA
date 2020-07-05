// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;
// using Numer::integrate;
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using std::log;
using std::exp;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

class comprisk: public Func
{
    private:
        MatrixXd lambda;
        VectorXd beta;
        VectorXd theta;
        VectorXd gamma;
        MatrixXd z;
        VectorXd w;
        VectorXd sj;
        VectorXd H;
        bool T_F; // T_F = 1: incidence rate for true response
        int G;
    public:
        comprisk(Rcpp::List list_, bool T_F_) :
            lambda( Rcpp::as<MapMatd>(list_["lambda"]) ),
            beta( Rcpp::as<MapVecd>(list_["beta"]) ),
            theta( Rcpp::as<MapVecd>(list_["theta"]) ),
            gamma( Rcpp::as<MapVecd>(list_["gamma"]) ),
            z( Rcpp::as<MapMatd>(list_["z"]) ),
            w( Rcpp::as<MapVecd>(list_["w"]) ),
            sj( Rcpp::as<MapVecd>(list_["sj"]) ),
            H( Rcpp::as<MapVecd>(list_["H"]) ),
            T_F( T_F_ ) {
            G = H.size();
        }

        double operator()(const double& t) const {
            VectorXd cum_lambda(2);
            cum_lambda.setZero();

            double res = 1.0;

            size_t g = 0;
            while (t > sj(g+1) && g < (G - 1) )  {
                cum_lambda += H(g) * lambda.col(g);
                g++;
            }
            for (size_t c=0; c<2; c++) {
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

        double hazard(const double &t, const bool cause) {
            double res;
            size_t g = 0;
            while (t > sj(g) && g < (G - 1) ) {
                g++;
            }
            if (cause) {
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
Rcpp::List cumcifun(Rcpp::List &param, bool T_F, double lower, double upper)
{
    comprisk f(param, T_F);
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
VectorXd cumcicurve(Rcpp::List &param, bool T_F, double lower, double upper, int npar)
{
    VectorXd res(npar);
    double delta = (upper - lower) / (double) npar;
    comprisk f(param, T_F);
    double err_est;
    int err_code;
    for (size_t n=0; n<npar; n++) {
        res(n) = integrate(f, lower, lower + (n + 1) * delta, err_est, err_code);
    }
    return res;
}

// [[Rcpp::export]]
VectorXd eval_incirate(Rcpp::List &param, bool T_F, VectorXd t)
{
    comprisk f(param, T_F);
    VectorXd res(t.size());
    for (size_t i=0; i<t.size(); i++) {
        res(i) = f(t(i));
    }
    return res;
}

// [[Rcpp::export]]
VectorXd eval_accuracy(Rcpp::List &param, VectorXd t)
{
    comprisk f(param, 1);
    VectorXd res(t.size());
    for (size_t i=0; i<t.size(); i++) {
        res(i) = f.accuracy(t(i));
    }
    return res;
}
