#include <Rcpp.h>

using namespace Rcpp;

// typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
// typedef Eigen::Map<Eigen::VectorXd> MapVecd;

NumericMatrix reshape(NumericMatrix::Row &row, int nrow, int ncol) {
  NumericMatrix res(nrow, ncol);
  for (int rr = 0; rr < nrow; rr++) {
    for (int cc = 0; cc < ncol; cc++) {
      res(rr, cc) = row[rr * ncol + cc];
    }
  }
  return res;
}

double distance(NumericMatrix::Row z_k, NumericMatrix::Row w_i) {
  double dist = 0;
    dist = sqrt(sum(pow((z_k - w_i), 2)));
  return dist;
}

class samples {
private:
  int I;
  int N;
  int G;
  // NumericVector sj;
  NumericMatrix H;
  NumericVector len;
  IntegerMatrix seg;
  IntegerMatrix Y;
  NumericMatrix lambda;
  NumericMatrix cum_lambda;
  NumericMatrix beta;
  NumericMatrix theta;
  NumericVector gamma;
  NumericMatrix z;
  NumericMatrix w;

public:
  samples(NumericMatrix::Row lambda_, NumericMatrix::Row theta_,
          NumericMatrix::Row z_, NumericMatrix::Row w_,
          NumericMatrix::Row gamma_, Rcpp::List param_)
      : I(as<int>(param_["I"])), N(as<int>(param_["N"])),
        G(as<int>(param_["G"])), // sj(as<NumericVector>(param_["sj"])),
        H(as<NumericMatrix>(param_["H"])),
        len(as<IntegerVector>(param_["len"])),
        seg(as<IntegerMatrix>(param_["seg"])), Y(as<IntegerMatrix>(param_["Y"])) {

    gamma[1] = gamma_[1];
    gamma[0] = gamma_[0];

    lambda = reshape(lambda_, 2 * I, G);
    theta = reshape(theta_, N, 2);
    z = reshape(z_, 2 * N, 2);
    w = reshape(w_, 2 * I, 2);

    // // Set number of rows and columns to attribute dim
    // lambda.attr("dim") = Dimension(G, 2 * I);
    // theta.attr("dim") = Dimension(2, N);
    // z.attr("dim") = Dimension(2, 2 * N);
    // w.attr("dim") = Dimension(2, 2 * I);

    // Converting to Rcpp Matrix type
    //lambda = transpose(lambda);
    //theta = transpose(theta);
    //z = transpose(z);
    //w = transpose(w);

    cum_lambda(2 * I, N); // matrx(0, 2 * I, N)

    for (int c = 0; c < 2; c++) {
      for (int i = 0; i < I; i++) {
        for (int k = 0; k < N; k++) {
          // if (NA(i, k) == 1) {
          for (int g = 0; g < seg(i, k); g++) {
            cum_lambda(i, k) += len(g) * lambda(i, g);
          }
          cum_lambda(i, k) += H(i, k) * lambda(i, seg(i, k));
          // }
        }
      }
    }
  }

  double loglike() {
    double running_total = 0;
    // std::cout << "Calculating the log-acceptance ratio...\n";
    for (int c = 0; c < 2; c++) {
      for (int i = 0; i < I; i++) {
        for (int k = 0; k < N; k++) {

          // if (NA(i, k) == 1) {
          // for (int g = 0; g < seg(i,k); g++) {
          //     cum_lambda(i,k) += len(g) * lambda(i,g);
          // }
          // cum_lambda(i,k) += H(i,k) * lambda(i,seg(i,k));

          running_total -=
              cum_lambda(c * I + i, k) *
              exp(beta(i, c) + theta(k, c) -
                  gamma(c) * distance(z.row(c * N + k), w.row(c * I + i)));

          if (Y(i, k) == c) {
            running_total +=
                log(lambda(c * I + i, seg(i, k))) + beta(i, c) + theta(k, c) -
                gamma(c) * distance(z.row(c * N + k), w.row(c * I + i));
          }
        }
      }
    }
    return running_total;
  }
};

// [[Rcpp::export]]
NumericVector get_loglike(NumericMatrix lambda_, NumericMatrix theta_,
                          NumericMatrix z_, NumericMatrix w_,
                          NumericMatrix gamma_, List &param_) {

  int num_iter = lambda_.nrow();
  NumericVector res(num_iter);

  for (int nn = 0; nn < num_iter; nn++) {
    samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
                     gamma_.row(nn), param_);
    res[nn] = sample_l.loglike();
  }
  return res;
}
