// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
#include <RcppParallel.h>
#include <RcppEigen.h>

using namespace Numer;
  
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::VectorXi> MapVeci;

class comprisk: public Func
{
    private:
        Eigen::MatrixXd lambda;
        Eigen::VectorXd beta;
        Eigen::VectorXd theta;
        Eigen::VectorXd gamma;
        Eigen::MatrixXd z;
        Eigen::MatrixXd w;
        Eigen::VectorXd sj;
        Eigen::VectorXd H;
        bool T_F; // T_F = 1: incidence rate for true response
        int G;
    public:
        comprisk(Rcpp::List list_, bool T_F_) :
            lambda( Rcpp::as<MapMatd>(list_["lambda"]) ),
            beta( Rcpp::as<MapVecd>(list_["beta"]) ),
            theta( Rcpp::as<MapVecd>(list_["theta"]) ),
            gamma( Rcpp::as<MapVecd>(list_["gamma"]) ),
            z( Rcpp::as<MapMatd>(list_["z"]) ),
            w( Rcpp::as<MapMatd>(list_["w"]) ),
            sj( Rcpp::as<MapVecd>(list_["sj"]) ),
            H( Rcpp::as<MapVecd>(list_["H"]) ),
            T_F( T_F_ ) {
            G = H.size();
        }

        double operator()(const double& t) const {
            Eigen::VectorXd cum_lambda(2);
            cum_lambda.setZero();

            double res = 1.0;

            int g = 0;
            while (t > sj(g+1) && g < (G - 1) )  {
                cum_lambda += H(g) * lambda.col(g);
                g++;
            }
            for (int c=0; c<2; c++) {
              res *=
                  std::exp(-1.0 * (cum_lambda(c) + (t - sj(g)) * lambda(c, g)) *
                           std::exp(beta(c) + theta(c) -
                                    gamma(c) * (z.row(c) - w.row(c)).norm()));
            }
            if (T_F) {
                res *= lambda(1,g) * std::exp( beta(1) + theta(1) - gamma(1) * (z.row(1)-w.row(1)).norm() );
            }
            else {
                res *= lambda(0,g) * std::exp( beta(0) + theta(0) - gamma(0) * (z.row(0)-w.row(0)).norm() );
            }
            return res;
        }

        double hazard(const double &t, const bool cause) {
            double res;
            int g = 0;
            while (t > sj(g) && g < (G - 1) ) {
                g++;
            }
            if (cause) {
                res = lambda(1,g) * std::exp( beta(1) + theta(1) - gamma(1) * (z.row(1)-w.row(1)).norm() );
            }
            else {
                res = lambda(0,g) * std::exp( beta(0) + theta(0) - gamma(0) * (z.row(0)-w.row(0)).norm() );
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
Eigen::VectorXd cumcicurve(Rcpp::List &param, bool T_F, double lower, double upper, int npar)
{
    Eigen::VectorXd res(npar);
    double delta = (upper - lower) / (double) npar;
    comprisk f(param, T_F);
    double err_est;
    int err_code;
    for (int n=0; n<npar; n++) {
        res(n) = integrate(f, lower, lower + (n + 1) * delta, err_est, err_code);
    }
    return res;
}

// [[Rcpp::export]]
Eigen::VectorXd eval_incirate(Rcpp::List &param, bool T_F, Eigen::VectorXd t)
{
    comprisk f(param, T_F);
    Eigen::VectorXd res(t.size());
    for (int i=0; i<t.size(); i++) {
        res(i) = f(t(i));
    }
    return res;
}

// [[Rcpp::export]]
Eigen::VectorXd eval_accuracy(Rcpp::List &param, Eigen::VectorXd t)
{
    comprisk f(param, 1);
    Eigen::VectorXd res(t.size());
    for (int i=0; i<t.size(); i++) {
        res(i) = f.accuracy(t(i));
    }
    return res;
}

// [[Rcpp::export]]
double eval_log_like(Rcpp::List &param, bool T_F, double t)
{
    comprisk f(param, T_F);
    double res;
    res = std::log(f(t));
    return res;
}
