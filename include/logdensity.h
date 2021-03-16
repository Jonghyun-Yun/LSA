#ifndef __LOGDENS_INCLUDED__
#define __LOGDENS_INCLUDED__

#include <stan/math.hpp>
#include <vector>

using namespace std;

class LogDensity {
  // Abstract base class for Log Densities
public:
  virtual vector<double> h(vector<double> x) = 0;
  virtual vector<double> h_prime(vector<double> x) = 0;
  virtual double h(double x) = 0;
  virtual double h_prime(double x) = 0;
};

class PoissonLogNormal : public LogDensity {
  // Log Density of poisson-log-normal posterior
  // log(pois(y; exp(z))) + log(norm(z; mu, sigma^2))
public:
  double y;
  double mu;
  double sigma;

  PoissonLogNormal(double y, double mu, double sigma);

  vector<double> clean(vector<double> x);

  vector<double> h(vector<double> x);
  vector<double> h_prime(vector<double> x);
  double h(double x);
  double h_prime(double x);
};

class gamma_full_conditional : public LogDensity {
    // Log Density of full-conditional for gamma
  public:
    int c;
    Eigen::VectorXd gamma;
    Eigen::VectorXd mu_gamma;
    Eigen::VectorXd sigma_gamma;
    Eigen::VectorXd jump_gamma;
    Eigen::MatrixXd cum_lambda;
    Eigen::MatrixXd beta;
    Eigen::MatrixXd theta;
    Eigen::MatrixXd z;
    Eigen::MatrixXd w;
    int I;
    int N;
    Eigen::MatrixXi NA;
    Eigen::MatrixXi Y;
    boost::ecuyer1988 rng;

    gamma_full_conditional(const int ONE_FREE_GAMMA, Eigen::VectorXd &gamma, Eigen::VectorXd &acc_gamma,
                           const Eigen::VectorXd &mu_gamma, const Eigen::VectorXd &sigma_gamma,
                           const Eigen::VectorXd &jump_gamma, const Eigen::MatrixXd &cum_lambda,
                           const Eigen::MatrixXd &beta, const Eigen::MatrixXd &theta,
                           const Eigen::MatrixXd &z, const Eigen::MatrixXd &w, const int &I,
                           const int &N, const Eigen::MatrixXi &NA,
                           const Eigen::MatrixXi &Y, boost::ecuyer1988 &rng);

    // vector<double> clean(vector<double> x);

    vector<double> h(vector<double> x);
    vector<double> h_prime(vector<double> x);
    double h(double x);
    double h_prime(double x);
};

#endif
