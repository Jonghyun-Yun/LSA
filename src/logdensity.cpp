#include "logdensity.h"
#include <cmath>
#include <stan/math.hpp>
#include "tbb/blocked_range2d.h"

using namespace std;

// PoissonLogNormal implementation
PoissonLogNormal::PoissonLogNormal(double y, double mu, double sigma) :
  y(y), mu(mu), sigma(sigma) {}

vector<double> PoissonLogNormal::clean(vector<double> x) {
  // ensure x is below 700 (max value for which we can compute exp)
  for (int i = 0; i < x.size(); ++i) {
    if (x.at(i) > 700) {
      x.at(i) = 700;
    }
  }
  return x;
}

vector<double> PoissonLogNormal::h(vector<double> x) {
  x = clean(x);
  vector<double> rv(x.size());
  for (int i = 0; i < x.size(); ++i) {
  rv.at(i) = h(x.at(i));
  }
  return rv;
}

vector<double> PoissonLogNormal::h_prime(vector<double> x) {
  x = clean(x);
  vector<double> rv(x.size());
  for (int i = 0; i < x.size(); ++i) {
    rv.at(i) = h_prime(x.at(i));
  }
  return rv;
}

double PoissonLogNormal::h(double x) {
  double rv = x*y - exp(x) - 0.5*pow(sigma, -2)*pow(x - mu, 2);
  return rv;
}

double PoissonLogNormal::h_prime(double x) {
  double rv =  y - exp(x) - pow(sigma, -2) * (x - mu);
  return rv;
}

// PoissonLogNormal implementation
gamma_full_conditional::gamma_full_conditional(const int ONE_FREE_GAMMA_, Eigen::VectorXd &gamma_, Eigen::VectorXd &acc_gamma_,
                                               const Eigen::VectorXd &mu_gamma_, const Eigen::VectorXd &sigma_gamma_,
                                               const Eigen::VectorXd &jump_gamma_, const Eigen::MatrixXd &cum_lambda_,
                                               const Eigen::MatrixXd &beta_, const Eigen::MatrixXd &theta_,
                                               const Eigen::MatrixXd &z_, const Eigen::MatrixXd &w_, const int &I_,
                                               const int &N_, const Eigen::MatrixXi &NA_,
                                               const Eigen::MatrixXi &Y_, boost::ecuyer1988 &rng_) :
  c(ONE_FREE_GAMMA_),
  gamma(gamma_),
  mu_gamma(mu_gamma_),
  sigma_gamma(sigma_gamma_),
  jump_gamma(jump_gamma_),
  cum_lambda(cum_lambda_),
  beta(beta_),
  theta(theta_),
  z(z_),
  w(w_),
  I(I_),
  N(N_),
  NA(NA_),
  Y(Y_),
  rng(rng_) {}

double gamma_full_conditional::h(double gamma_s) {
  double rv = stan::math::lognormal_lpdf(gamma_s, mu_gamma(c), sigma_gamma(c));
  rv += tbb::parallel_reduce(
    tbb::blocked_range2d<int>(0, I, 0, N), 0.0,
    [&](tbb::blocked_range2d<int> r, double running_total) {
      for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
        for (int k = r.cols().begin(); k < r.cols().end(); ++k) {
          if (NA(i, k) == 1) {

            running_total -=
              cum_lambda(c * I + i, k) *
              exp(beta(i, c) + theta(k, c) -
                  gamma_s * stan::math::distance(z.row(c * N + k),
                                           w.row(c * I + i)));

            if (Y(i, k) == c) {
              running_total -=
                gamma_s * stan::math::distance(z.row(c * N + k), w.row(c * I + i));
            }
          }
        }
      }
      return running_total;
    },
    std::plus<double>());

  return rv;
}

double gamma_full_conditional::h_prime(double gamma_s) {
  double rv = - 1.0 * (1.0 + (log(gamma_s) - mu_gamma(c)) / pow(sigma_gamma(c), 2)) / gamma_s; // lognormal_lpdf derivative
  rv += tbb::parallel_reduce(
    tbb::blocked_range2d<int>(0, I, 0, N), 0.0,
    [&](tbb::blocked_range2d<int> r, double running_total) {
      for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
        for (int k = r.cols().begin(); k < r.cols().end(); ++k) {
          if (NA(i, k) == 1) {

            running_total +=
              cum_lambda(c * I + i, k) *
              exp(beta(i, c) + theta(k, c) -
                  gamma_s * stan::math::distance(z.row(c * N + k),
                                           w.row(c * I + i))) * stan::math::distance(z.row(c * N + k),
                                           w.row(c * I + i));

            if (Y(i, k) == c) {
              running_total -= stan::math::distance(z.row(c * N + k), w.row(c * I + i));
            }
          }
        }
      }
      return running_total;
    },
    std::plus<double>());

  return rv;
}

vector<double> gamma_full_conditional::h(vector<double> gamma_s) {
  vector<double> rv(gamma_s.size());
  for (int i = 0; i < gamma_s.size(); ++i) {
    rv.at(i) = h(gamma_s.at(i));
  }
  return rv;
}

vector<double> gamma_full_conditional::h_prime(vector<double> gamma_s) {
  vector<double> rv(gamma_s.size());
  for (int i = 0; i < gamma_s.size(); ++i) {
    rv.at(i) = h_prime(gamma_s.at(i));
  }
  return rv;
}
