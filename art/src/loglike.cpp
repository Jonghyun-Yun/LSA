// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <algorithm>

using MapMatd = Eigen::Map<Eigen::MatrixXd>;
using MapVecd = Eigen::Map<Eigen::VectorXd>;
// typedef Eigen::Map<Eigen::MatrixXi> MapMati; // cannot map integermatrix??
using MapVeci = Eigen::Map<Eigen::VectorXi>;

using Rcpp::as;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

// find i such that x is in i-th interval (for i=1,...,(# breaks - 1))
// x below the lower bound gives i = 0
// x above the upper bound gives i = # of breaks
int findInterval(const double x, const VectorXd &breaks) noexcept {
  int out = 0;
  for (int i = 0; i < breaks.size(); i++) {
    if (x <= breaks(i))
      break;
    ++out;
  }
  return out;
}

template <typename T>
VectorXi findInterval(const T &x, const VectorXd &breaks) noexcept {
  VectorXi out(x.size());

  for (int i = 0; i < x.size(); i++) {
    out(i) = findInterval(x(i), breaks);
  }

  return out;
}

// first elem 0
VectorXd my_cumsum(const VectorXd &x) {
  // initialize the result vector
  VectorXd res(x.size() + 1);
  res(0) = 0;
  // std::partial_sum(x.begin(), x.end(), res.begin());
  for (int i = 1; i <= x.size(); i++) {
    res(i) = res(i - 1) + x(i - 1);
  }
  // std::partial_sum(x.begin(), x.end(), res.begin());
  return res;
}

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

// 2D l2-norm
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
  VectorXd sj;
  MatrixXd H;
  VectorXi len;
  MatrixXi seg;
  MatrixXi Y;
  MatrixXd lambda;
  MatrixXd cum_lambda;
  MatrixXd beta;
  MatrixXd theta;
  Eigen::Vector2d gamma;
  MatrixXd z;
  MatrixXd w;

public:
  samples(NumericMatrix::Row lambda_, NumericMatrix::Row theta_,
          NumericMatrix::Row z_, NumericMatrix::Row w_,
          NumericMatrix::Row gamma_, Rcpp::List param_)
      : //   : I(as<int>(param_["I"])), N(as<int>(param_["N"])),
        //     G(as<int>(param_["G"])) {
        sj(as<MapVecd>(param_["sj"])), H(as<MapMatd>(param_["H"])),
        len(as<MapVeci>(param_["len"])) {

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

  void set_seg(int i, const VectorXi &seg_i) {
    try {
      if (seg_i.size() != N) throw 0;
    }
    catch (int n) {std::cout << "size mismatched" << std::endl;}
    seg.row(i) = seg_i;
  }

  void eval_cum_lambda() noexcept {
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
    return theta(k, c) -
           gamma(c) * (z.row(c * N + k) - w.row(c * I + i)).norm();
  }

  // log(lambda) + log_rr
  double log_hazard(int c, int i, int k) {
    return std::log(lambda(c * I + i, seg(i, k))) + log_rr(c, i, k);
  }

  double hazard(int c, int i, int k) { return std::exp(log_hazard(c, i, k)); }

  double eval_acc(int i, int k) {
    return hazard(1, i, k) / (hazard(1, i, k) + hazard(0, i, k));
  }

  // void print_param(int i, int k) {
  //   for (int c = 0; c < 2; c++) {
  //   std::cout << seg(i,k) << "," << lambda(c * I + i,seg(i,k)) << "," <<
  //   theta(k,c) << "," << gamma(c) << "," << (z.row(c * N + k) - w.row(c * I +
  //   i)).norm() << "," << hazard(c,i,k) << std::endl; std::cout <<
  //   "====================================\n";
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

  double loglike() {
    eval_cum_lambda();
    double running_total = 0;
    //   // std::cout << "Calculating the log-acceptance ratio...\n";
    for (int c = 0; c < 2; c++) {
      for (int i = 0; i < I; i++) {
        for (int k = 0; k < N; k++) {

          running_total -= cum_lambda(c * I + i, k) * std::exp(log_rr(c, i, k));

          if (Y(i, k) == c) {
            running_total += log_hazard(c, i, k);
          }
        }
      }
    }
    return running_total;
  }

  double gen_time(int i, int k) {
    double vtime;
    VectorXd haz = VectorXd::Zero(G);
    VectorXd hazlen = VectorXd::Zero(G);

    for (int g = 0; g < G; g++) {
      for (int c = 0; c < 2; c++) {
        haz(g) += lambda(c * I + i, g) * std::exp(log_rr(c, i, k));
      }
      hazlen(g) = haz(g) * len(g);
    }

    VectorXd cumhaz = my_cumsum(hazlen);

    // double logS = -1.0;
    double logS = -1.0 * R::rexp(1.0);
    // int mj = G - 1;
    int ss = findInterval(-1.0 * logS, cumhaz) - 1;

    // try {
    //   if ((haz.array() < 0).any())
    //     throw 0;
    // } catch (int n) {
    //   std::cout << "Caught " << n << std::endl;
    //       std::cout << "i, k: " << i << "," << k << std::endl;
    //       for (int c = 0; c < 2; c++) {
    // std::cout <<
    //       theta(k, c) << "," << gamma(c) << "," << (z.row(c * N + k) -
    //       w.row(c * I + i)).norm() << std::endl;

    // std::cout << z.row(c * N + k) << "," << w.row(c * I + i)<< std::endl;
    //         for (int g = 0; g < G; g++) {

    //           std::cout << "c, g: " << c << "," << g << std::endl;
    //           std::cout << "lambda(c * I + i, g)" << lambda(c * I + i, g)
    //                     << std::endl;
    //           std::cout << "rr" << std::exp(log_rr(c, i, k)) << std::endl;
    //           std::cout << "len" << len(g) << std::endl;
    //         }
    //       }

    // std::cout << "cumhaz" << std::endl;
    // std::cout << cumhaz << std::endl;
    // std::cout << "-1.0 * logS" << std::endl;
    // std::cout << -1.0 * logS << std::endl;
    // std::cout << "ss" << std::endl;
    // std::cout << ss << std::endl;

    // }

    if (ss < G) {
      vtime = sj(ss) - (logS + cumhaz(ss)) / haz(ss);
    } else
      vtime = sj(ss);
    return vtime;
  }

  VectorXd gen_vtime(int item) noexcept {
    VectorXd res(N);
    for (int k = 0; k < N; k++) {
      res(k) = gen_time(item, k);
    }
    return res;
  }

  double logprior() { return 0; }
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

// TODO: create member functions
// [[Rcpp::export]]
Eigen::MatrixXd rcpp_gen_surv_pp(NumericMatrix lambda_, NumericMatrix theta_,
                                 NumericMatrix z_, NumericMatrix w_,
                                 NumericMatrix gamma_, List &param_, int item) {

  int num_iter = lambda_.nrow();
  int N = theta_.ncol() / 2;
  int I = w_.nrow() / 2;

  try {
    if (item < 1 || item > I)
      throw 0;
  } catch (int n) {
    std::cout << "item must be between 1 and I." << std::endl;
  }
  item -= 1; // R to C indexing

  Eigen::MatrixXd res(num_iter, N);

  // std::cout << "res: " << res.at(1) << std::endl;

  for (int nn = 0; nn < num_iter; nn++) {
    samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
                     gamma_.row(nn), param_);
    res.row(nn) = sample_l.acc_item(item);
    Rcpp::checkUserInterrupt();
  }
  // int nn = 0;
  // samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
  //                  gamma_.row(nn), param_);
  // sample_l.print_param(39, 0);
  return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_gen_surv_time(NumericMatrix lambda_, NumericMatrix theta_,
                                   NumericMatrix z_, NumericMatrix w_,
                                   NumericMatrix gamma_, List &param_,
                                   int item) {

  int num_iter = lambda_.nrow();
  int N = theta_.ncol() / 2;
  int I = w_.nrow() / 2;

  try {
    if (item < 1 || item > I)
      throw 0;
  } catch (int n) {
    std::cout << "item must be between 1 and I." << std::endl;
  }
  item -= 1; // R to C indexing

  Eigen::MatrixXd res(num_iter, N);

  // std::cout << "res: " << res.at(1) << std::endl;

  for (int nn = 0; nn < num_iter; nn++) {
    samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
                     gamma_.row(nn), param_);
    res.row(nn) = sample_l.gen_vtime(item);
    Rcpp::checkUserInterrupt();
  }
  // int nn = 0;
  // samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
  //                  gamma_.row(nn), param_);
  // sample_l.print_param(39, 0);
  return res;
}

// calculate segment for new time values
template <typename T>
VectorXi time_seg(const T &x, const VectorXd &sj) noexcept {
  VectorXi out(x.size());

  for (int i = 0; i < x.size(); i++) {
    out(i) = findInterval(x(i), sj) - 1;
    if (out(i) == (sj.size() - 1))
      out(i) -= 1; // out of range
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List rcpp_gen_surv(NumericMatrix lambda_, NumericMatrix theta_,
                         NumericMatrix z_, NumericMatrix w_,
                         NumericMatrix gamma_, List &param_, int item) {

  int num_iter = lambda_.nrow();
  int N = theta_.ncol() / 2;
  int I = w_.nrow() / 2;
  // int I = w_.nrow() / 2;
  VectorXd sj(as<MapVecd>(param_["sj"]));

  try {
    if (item < 1 || item > I)
      throw 0;
  } catch (int n) {
    std::cout << "item must be between 1 and I." << std::endl;
  }
  item -= 1; // R to C indexing

  // if (item < 0 || item > (I-1)) return 0;

  Eigen::MatrixXd res_t(num_iter, N);
  Eigen::MatrixXd res_p(num_iter, N);

  // std::cout << "res: " << res.at(1) << std::endl;
  for (int nn = 0; nn < num_iter; nn++) {
    // for (int nn = 0; nn < 1; nn++) {
    samples sample_l(lambda_.row(nn), theta_.row(nn), z_.row(nn), w_.row(nn),
                     gamma_.row(nn), param_);
    res_t.row(nn) = sample_l.gen_vtime(item); // indexing adjust
    sample_l.set_seg(item, time_seg(res_t.row(nn), sj));
    res_p.row(nn) = sample_l.acc_item(item); // indexing adjust
    Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(Rcpp::Named("time") = res_t,
                            Rcpp::Named("pp") = res_p);
}
