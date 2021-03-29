// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
// typedef Eigen::Map<Eigen::MatrixXi> MapMati; // cannot map integermatrix??
typedef Eigen::Map<Eigen::VectorXi> MapVeci;

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::as;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;
using Eigen::VectorXi;

// maping rowvector to matrix by row-major order
MatrixXd reshape_rowmajor(NumericMatrix::Row &row, int nrow, int ncol) {
  MatrixXd res(nrow, ncol);
  for (int rr = 0; rr < nrow; rr++) {
    for (int cc = 0; cc < ncol; cc++) {
      res(rr, cc) = row[rr * ncol + cc];
    }
  }
  return res;
}

double distance(NumericMatrix::Row z_k, NumericMatrix::Row w_i) {
  double dist = 0;
  for (int d = 0; d < 2; d++) {
    dist += std::pow(z_k[d] - w_i[d], 2);
  }
  return std::sqrt(dist);
}

class samples {
private:
  int I;
  int N;
  int G;
  // NumericVector sj;
  MatrixXd H;
  VectorXi len;
  MatrixXi seg;
  MatrixXi Y;

public:
  MatrixXd lambda;
  MatrixXd cum_lambda;
  MatrixXd beta;
  MatrixXd theta;
  Eigen::Vector2d gamma;
  MatrixXd z;
  MatrixXd w;

  samples(NumericMatrix::Row lambda_, NumericMatrix::Row theta_,
          NumericMatrix::Row z_, NumericMatrix::Row w_,
          NumericMatrix::Row gamma_, Rcpp::List param_):
    //   : I(as<int>(param_["I"])), N(as<int>(param_["N"])),
    //     G(as<int>(param_["G"])) {
    //   // sj(as<NumericVector>(param_["sj"])),
        H(as<MapMatd>(param_["H"])),
        len(as<MapVeci>(param_["len"]))
    {

    MatrixXd tmp_seg = as<MapMatd>(param_["seg"]);
    MatrixXd tmp_Y = as<MapMatd>(param_["Y"]);

    seg = tmp_seg.cast<int>();
    Y = tmp_Y.cast<int>();

    I = H.rows();
    N = H.cols();
    G = len.size();

    // std::cout << "G" << std::endl;
    // std::cout << G << std::endl;

    // VectorXd gamma(2);
    gamma(0) = gamma_[0];
    gamma(1) = gamma_[1];

    lambda = reshape_rowmajor(lambda_, 2 * I, G);
    theta = reshape_rowmajor(theta_, N, 2);
    z = reshape_rowmajor(z_, 2 * N, 2);
    w = reshape_rowmajor(w_, 2 * I, 2);
    }

    void eval_cum_lambda() {
      cum_lambda.setZero(2 * I, N); // matrx(0, 2 * I, N)

      for (int c = 0; c < 2; c++) {
        for (int i = 0; i < I; i++) {
          for (int k = 0; k < N; k++) {
            // if (NA(i, k) == 1) {
            for (int g = 0; g < seg(i, k); g++) {
              // std::cout << seg(i,k) << std::endl;
              cum_lambda(c * I + i, k) += len(g) * lambda(c * I + i, g);
            }
            cum_lambda(c * I + i, k) += H(i, k) * lambda(c * I + i, seg(i, k));
            // }
          }
        }
      }
    }


    double log_rr(int c, int i, int k) {
      double lr;
      lr = theta(k, c) -
                  gamma(c) * (z.row(c * N + k) - w.row(c * I + i)).norm();
      return lr;
    }

    double log_hazard(int c, int i, int k) {
      double lh;
      lh = std::log(lambda(c * I + i, seg(i, k))) + log_rr(c,i,k);
      return lh;
    }

    double hazard(int c, int i, int k) {
      double h;
      h = std::exp(log_hazard(c,i,k));
      return h;
    }

    double eval_acc(int i, int k) {
      double acc;
      acc = hazard(1,i,k) / (hazard(1,i,k) + hazard(0,i,k));
      return acc;
    }

    // void print_param(int i, int k) {
    //   for (int c = 0; c < 2; c++) {
    //   std::cout << seg(i,k) << "," << lambda(c * I + i,seg(i,k)) << "," << theta(k,c) << "," << gamma(c) << "," << (z.row(c * N + k) - w.row(c * I + i)).norm() << "," << hazard(c,i,k) << std::endl;
    //   std::cout << "====================================\n";
    //   }
    //   std::cout << eval_acc(i,k) << std::endl;;
    // }

    VectorXd acc_item(int item) {
      VectorXd acc(N);
      for (int k = 0; k < N; k++) {
      acc(k) = eval_acc(item, k);
      }
      return acc;
    }

    void new_seg(int i, VectorXi seg_i) {
     seg.col(i) = seg_i;
    }

  double loglike() {
    eval_cum_lambda();
    double running_total = 0;
  //   // std::cout << "Calculating the log-acceptance ratio...\n";
    for (int c = 0; c < 2; c++) {
      for (int i = 0; i < I; i++) {
        for (int k = 0; k < N; k++) {

          running_total -=
              cum_lambda(c * I + i, k) *
              std::exp(log_rr(c,i,k));

          if (Y(i, k) == c) {
            running_total += log_hazard(c,i,k);
          }
        }
      }
    }
    return running_total;
  }

    double logprior() {
      return 0;
    }
};

// [[Rcpp::export]]
NumericVector get_loglike(NumericMatrix lambda_, NumericMatrix theta_,
                          NumericMatrix z_, NumericMatrix w_,
                          NumericMatrix gamma_, List &param_) {

  int num_iter = lambda_.nrow();
  NumericVector res(num_iter);

  // std::cout << "res: " << res.at(1) << std::endl;

  for (int nn = 0; nn < num_iter; nn++) {
    samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
                     gamma_.row(nn), param_);
    res[nn] = sample_l.loglike();
  }
  return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_gen_surv_pp(NumericMatrix lambda_, NumericMatrix theta_,
                            NumericMatrix z_, NumericMatrix w_,
                            NumericMatrix gamma_, List &param_, int item) {

  int num_iter = lambda_.nrow();
  int N = theta_.ncol() / 2;
  // int I = w_.nrow() / 2;

  // if (item < 0 || item > (I-1)) return 0;

  Eigen::MatrixXd res(num_iter, N);

  // std::cout << "res: " << res.at(1) << std::endl;

  for (int nn = 0; nn < num_iter; nn++) {
    samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
                     gamma_.row(nn), param_);
    res.row(nn) = sample_l.acc_item(item-1); // indexing adjust
    Rcpp::checkUserInterrupt();
  }
  // int nn = 0;
  // samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
  //                  gamma_.row(nn), param_);
  // sample_l.print_param(39, 0);
  return res;
}
