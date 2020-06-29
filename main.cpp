#include <stan/math.hpp>
#include <iostream>
#include "my_header.h"
// basic file operations
#include <fstream>
// #include <cmath>
#include <cstdlib>  // to process main(arguemnts)
#include <ctime>    // std::time
#include <typeinfo> // typeid(a).name()

// #include "is_empty.h"
#include "create_rng.hpp"
#include "readCSV.h"
#include "update_lambda.h"
#include "update_theta.h"
#include "update_beta.h"
#include "update_z.h"
#include "update_w.h"
#include "update_gamma.h"
#include "fun_lp.h"

// see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
const static Eigen::IOFormat CSVFormatN(Eigen::StreamPrecision,
                                        Eigen::DontAlignCols, ", ", "\n");

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, ", ", ", ");

const static Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
                                          Eigen::DontAlignCols, ", ", ", ", "",
                                          "", "", "\n");

inline double update_sigma(const Eigen::MatrixXd &theta, const int &N, boost::ecuyer1988 &rng) {
  return stan::math::sqrt(stan::math::inv_gamma_rng(1.0 + 2.0 * N, 1.0 + theta.sum() / 2.0, rng));
}

int main(int argc, const char *argv[]) {

  int num_samples = atoi(argv[1]);
  int num_warmup = atoi(argv[2]);
  int thin = atoi(argv[3]);
  // double my_eps = atof(argv[4]);
  // int my_L = atoi(argv[5]);
  // double min_E = atof(argv[6]);
  // double max_E = atof(argv[7]);

  double lp_;

  int num_iter = num_samples + num_warmup;
  int num_print = num_iter / 100;
  int chain_id = 1; // chain ID

  std::ofstream osummary, osample;
  std::stringstream fstart, fsummary, fsample, fwarmup;

  fstart << "./output/start.csv";
  fwarmup << "./output/warmup.csv";
  fsummary << "./ouput/summary.csv";
  fsample << "./output/sample.csv";

  // use current time as seed for random generator
  int rseed = (unsigned int)time(0) / 2;
  boost::ecuyer1988 rng = create_rng(rseed, chain_id);
  (void)rng; // suppress unused var warning

  std::cout << std::fixed << std::setprecision(1);
  osummary.open(fstart.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }
  // if (is_empty(osummary)) inp << ".chain, " << "seed, " << "iter, " << "warm,
  // " << "thin" << std::endl;
  osummary << chain_id << ", " << rseed << ", " << num_samples << ", " << num_warmup
           << ", " << thin << std::endl;
  osummary.close();

  std::cout << "Reading data sets...\n";
  Eigen::MatrixXd mvar(1, 4);
  mvar = readCSV("mvar.csv", 1, 4);
  const int I = mvar(0, 0);
  const int N = mvar(0, 1);
  const int C = mvar(0, 2);
  const int G = mvar(0, 3);

  Eigen::MatrixXd tmp_mseg(I, N);
  Eigen::MatrixXd tmp_mY(I, N);
  Eigen::VectorXd mlen(G);
  Eigen::MatrixXi mseg(I, N);
  Eigen::MatrixXd mH(I, N);
  Eigen::MatrixXd mt(I, N);
  Eigen::MatrixXi mY(I, N);

  mlen = readCSV("mlen.csv", G, 1);
  tmp_mseg = readCSV("mseg.csv", I, N);
  mH = readCSV("mh.csv", I, N);
  mt = readCSV("mt.csv", I, N);
  tmp_mY = readCSV("mi.csv", I, N);

  mseg = tmp_mseg.cast<int>();
  mY = tmp_mY.cast<int>();

  std::cout << "Reading hyperparameters...\n";
  Eigen::MatrixXd tmp_lambda(3*I, G);
  Eigen::MatrixXd tmp_beta(3*I, 2);
  Eigen::MatrixXd tmp_theta(3*N, 2);
  Eigen::VectorXd tmp_sigma(2);
  Eigen::MatrixXd tmp_gamma(3, 2);
  Eigen::MatrixXd tmp_z(3*N, 2);
  Eigen::MatrixXd tmp_w(3*I, 2);

  std::cout << "Reading hyperparameters for lambda...\n";
  tmp_lambda = readCSV("pj_lambda.csv", 3*I, G);
  std::cout << "Reading hyperparameters for beta...\n";
  tmp_beta = readCSV("pj_beta.csv", 3*I, 2);
  std::cout << "Reading hyperparameters for theta...\n";
  tmp_theta = readCSV("pj_theta.csv", 3*N, 2);
  std::cout << "Reading hyperparameters for sigma...\n";
  tmp_sigma = readCSV("pj_sigma.csv", 2, 1);
  std::cout << "Reading hyperparameters for gamma...\n";
  tmp_gamma = readCSV("pj_gamma.csv", 3, 2);
  std::cout << "Reading hyperparameters for z...\n";
  tmp_z = readCSV("pj_z.csv", 3*N, 2);
  std::cout << "Reading hyperparameters for w...\n";
  tmp_w = readCSV("pj_w.csv", 3*I, 2);

  std::cout << "Reading hyperparameters for lambda...\n";
  Eigen::MatrixXd a_lambda(I,G);
  Eigen::MatrixXd b_lambda(I,G);
  Eigen::MatrixXd jump_lambda(I,G);

  a_lambda = tmp_lambda.topRows(I);
  b_lambda = tmp_lambda.block(I,0,I,G);
  jump_lambda = tmp_lambda.bottomRows(I);

  std::cout << "Reading hyperparameters for beta...\n";
  Eigen::MatrixXd mu_beta(I,2);
  Eigen::MatrixXd sigma_beta(I,2);
  Eigen::MatrixXd jump_beta(I,2);

  mu_beta = tmp_beta.topRows(I);
  sigma_beta = tmp_beta.block(I,0,I,2);
  jump_beta = tmp_beta.bottomRows(I);

  std::cout << "Reading hyperparameters for theta...\n";
  Eigen::MatrixXd mu_theta(N,2);
  Eigen::MatrixXd sigma_theta(N,2);
  Eigen::MatrixXd jump_theta(N,2);

  mu_theta = tmp_theta.topRows(N);
  sigma_theta = tmp_theta.block(N,0,N,2);
  jump_theta = tmp_theta.bottomRows(N);
  
  std::cout << "Reading hyperparameters for sigma...\n";
  double a_sigma = tmp_sigma(0);
  double b_sigma = tmp_sigma(1);

  std::cout << "Reading hyperparameters for gamma...\n";
  Eigen::VectorXd mu_gamma(2);
  Eigen::VectorXd sigma_gamma(2);
  Eigen::VectorXd jump_gamma(2);

  // mu_gamma = tmp_gamma.topRows(1);
  // sigma_gamma = tmp_gamma.row(1);
  // jump_gamma = tmp_gamma.bottomRows(1);

  std::cout << "Reading hyperparameters for z...\n";
  Eigen::MatrixXd mu_z(N,2);
  Eigen::MatrixXd sigma_z(N,2);
  Eigen::MatrixXd jump_z(N,2);

  mu_z = tmp_z.topRows(N);
  sigma_z = tmp_z.block(N,0,N,2);
  jump_z = tmp_z.bottomRows(N);
  
  std::cout << "Reading hyperparameters for w...\n";
  Eigen::MatrixXd mu_w(I,2);
  Eigen::MatrixXd sigma_w(I,2);
  Eigen::MatrixXd jump_w(I,2);
  
  mu_w = tmp_w.topRows(I);
  sigma_w = tmp_w.block(I,0,I,2);
  jump_w = tmp_w.bottomRows(I);
  
  std::cout << "Initializing parameters...\n";
  // parameters
  Eigen::MatrixXd theta(N, 2);
  Eigen::MatrixXd beta(I, 2);
  double sigma;
  Eigen::MatrixXd z0(N, 2);
  Eigen::MatrixXd z1(N, 2);
  Eigen::MatrixXd w(I, 2);
  Eigen::VectorXd gamma(2);
  Eigen::MatrixXd lambda0(I, G);
  Eigen::MatrixXd lambda1(I, G);

  // Initialization
  theta.setZero();
  beta.setZero();
  sigma = 1;
  z0.setZero();
  z1.setZero();
  w.setZero();
  gamma.setOnes();
  lambda0.setOnes();
  lambda1.setOnes();

  // accept_stat
  Eigen::MatrixXd acc_theta(N, 2);
  Eigen::MatrixXd acc_beta(I, 2);
  double acc_sigma;
  Eigen::VectorXd acc_z0(N);
  Eigen::VectorXd acc_z1(N);
  Eigen::MatrixXd acc_w(I, 2);
  Eigen::VectorXd acc_gamma(2);
  Eigen::MatrixXd acc_lambda0(I, G);
  Eigen::MatrixXd acc_lambda1(I, G);

  acc_theta.setZero();
  acc_beta.setZero();
  acc_sigma = 1;
  acc_z0.setZero();
  acc_z1.setZero();
  acc_w.setZero();
  acc_gamma.setZero();
  acc_lambda0.setZero();
  acc_lambda1.setZero();

  std::cout << "Starting Sampling...\n";
  for (int ii = 1; ii <= num_warmup; ii++) {

    std::clock_t c_start = std::clock();

    // updating lambda...
    for (int i = 0; i < I; i++) {
      for (int g = 0; g < G; g++) {
        update_lambda(lambda0(i, g), acc_lambda0(i, g), g,
                      beta(i, 0), theta.col(0), gamma(0), z0, w.row(i),
                      N, mlen(g), mseg.row(i), mH.row(i), mY.row(i), 0, rng);

        update_lambda(lambda1(i, g), acc_lambda1(i, g), g,
                      beta(i, 1), theta.col(1), gamma(1), z1, w.row(i),
                      N, mlen(g), mseg.row(i), mH.row(i), mY.row(i), 1, rng);
      }
    }

    // updating theta...
    for (int k = 0; k < N; k++) {

    update_theta( theta(k,0), acc_theta(k,0), sigma,
                  lambda0, beta.col(0), gamma(0), z0.row(k), w,
                  I, mlen, mseg.col(k), mH.col(k), mY.col(k), 0, rng);

    update_theta( theta(k,1), acc_theta(k,1), sigma,
                  lambda1, beta.col(1), gamma(1), z1.row(k), w,
                  I, mlen, mseg.col(k), mH.col(k), mY.col(k), 1, rng);

    }

    // updating beta...
    for (int i = 0; i < I; i++) {

      update_beta(beta(i,0), acc_beta(i,0),
                  lambda0.row(i), theta.col(0), gamma(0), z0, w.row(i),
                  N, mlen, mseg.row(i), mH.row(i), mY.row(i), 0, rng);
      
      update_beta(beta(i,1), acc_beta(i,1),
                  lambda1.row(i), theta.col(1), gamma(1), z1, w.row(i),
                  N, mlen, mseg.row(i), mH.row(i), mY.row(i), 1, rng);

    }

    // updating z...
    for (int k = 0; k < N; k++) {

      z0.row(k) = update_z(z0.row(k), acc_z0(k),
                           lambda0, beta.col(0), theta(k,0), gamma(0), w,
                           I, mlen, mseg.col(k), mH.col(k), mY.col(k), 0, rng);

      z1.row(k) = update_z(z1.row(k), acc_z1(k),
                           lambda1, beta.col(1), theta(k,1), gamma(1), w,
                           I, mlen, mseg.col(k), mH.col(k), mY.col(k), 1, rng);

    }

    // updating w...
    for (int i = 0; i < I; i++) {

      w.row(i) = update_w(w.row(i), acc_w(i),
                          lambda0.row(i), lambda1.row(i), beta.row(i), theta, gamma, z0, z1,
                          N, G, mlen, mseg.row(i), mH.row(i), mY.row(i), rng);

    }

    // updating gamma...
    update_gamma(gamma, acc_gamma,
                 lambda0, lambda1, beta, theta, z0, z1, w,
                 I, N, G, mlen, mseg, mH, mY, rng); 

    // updating sigma 
    sigma = update_sigma(theta, N, rng);
  
    if (ii == 10) {
      std::clock_t c_end = std::clock();
      double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
      std::cout << "\n10 iterations would take " << time_elapsed_ms
                << " milliseconds.\n"
                << "Adjust your expectations accordingly!\n\n";
    }
    if ((ii % num_print == 0) && (ii > 10))
      std::cout << "Iteration: " << std::setw(7) << std::right << ii << " / "
                << std::setw(7) << std::right << num_iter << " ["
                << std::setw(3) << (int)(ii / (double)num_iter * 100) << "%]"
                << "  (Warmup)" << std::endl;
  }

  // writing (overall) acceptance rate... (last column for sigma's acceptance rate = 1)
  osummary.open(fwarmup.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }
  // if (is_empty(osummary)) osummary << ".chain, " << "accept_stat\n";
  osummary << chain_id << ", " << acc_lambda0.mean() / (double)num_warmup
           << ", " << acc_lambda1.mean() / (double)num_warmup
           << ", " << acc_beta.mean() / (double)num_warmup
           << ", " << acc_theta.mean() / (double)num_warmup
           << ", " << acc_z0.mean() / (double)num_warmup
           << ", " << acc_z1.mean() / (double)num_warmup
           << ", " << acc_w.mean() / (double)num_warmup
           << ", " << acc_gamma.mean() / (double)num_warmup
           << ", " << acc_sigma
           << std::endl;
  osummary.close();

  std::cout << std::fixed << std::setprecision(1);
  osample.open(fsample.str(), std::ios::app);
  if (!osample.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }

  // init: acceptance ratio
  acc_theta.setZero();
  acc_beta.setZero();
  acc_sigma = 1;
  acc_z0.setZero();
  acc_z1.setZero();
  acc_w.setZero();
  acc_gamma.setZero();
  acc_lambda0.setZero();
  acc_lambda1.setZero();

  for (int nn = 1; nn <= num_samples; nn++) {

    // lambda updating
    for (int i = 0; i < I; i++) {
      for (int g = 0; g < G; g++) {
        update_lambda(lambda0(i, g), acc_lambda0(i, g), g, beta(i, 0),
                      theta.col(0), gamma(0), z0, w.row(i),
                      N, mlen(g), mseg.row(i), mH.row(i), mY.row(i), 0, rng);

        update_lambda(lambda1(i, g), acc_lambda1(i, g), g,
                      beta(i, 1), theta.col(1), gamma(1), z1, w.row(i),
                      N, mlen(g), mseg.row(i), mH.row(i), mY.row(i), 1, rng);
      }
    }

    // theta updating

    for (int k = 0; k < N; k++) {

      update_theta( theta(k,0), acc_theta(k,0), sigma,
                    lambda0, beta.col(0), gamma(0), z0.row(k), w,
                    I, mlen, mseg.col(k), mH.col(k), mY.col(k), 0, rng);

      update_theta( theta(k,1), acc_theta(k,1), sigma,
                    lambda1, beta.col(1), gamma(1), z1.row(k), w,
                    I, mlen, mseg.col(k), mH.col(k), mY.col(k), 1, rng);

    }

    // updating beta...
    for (int i = 0; i < I; i++) {

      update_beta(beta(i,0), acc_beta(i,0),
                  lambda0.row(i), theta.col(0), gamma(0), z0, w.row(i),
                  N, mlen, mseg.row(i), mH.row(i), mY.row(i), 0, rng);
      
      update_beta(beta(i,1), acc_beta(i,1),
                  lambda1.row(i), theta.col(1), gamma(1), z1, w.row(i),
                  N, mlen, mseg.row(i), mH.row(i), mY.row(i), 1, rng);

    }

    // updating z...
    for (int k = 0; k < N; k++) {

      z0.row(k) = update_z(z0.row(k), acc_z0(k),
                           lambda0, beta.col(0), theta(k,0), gamma(0), w,
                           I, mlen, mseg.col(k), mH.col(k), mY.col(k), 0, rng);

      z1.row(k) = update_z(z1.row(k), acc_z1(k),
                           lambda1, beta.col(1), theta(k,1), gamma(1), w,
                           I, mlen, mseg.col(k), mH.col(k), mY.col(k), 1, rng);

    }
    
    if ((nn % num_print == 0) || nn ==  num_samples)
      std::cout << "Iteration: " << std::setw(7) << std::right << (nn + num_warmup) << " / "
                << std::setw(7) << std::right << num_iter << " ["
                << std::setw(3)
                << (int)((nn + num_warmup) / (double)num_iter * 100) << "%]"
                << "  (Sampling)" << std::endl;
    if (nn % thin == 0) {
// eval log_prob
      lp_ = fun_lp(a_lambda, b_lambda, mu_beta, sigma_beta, mu_theta, sigma_theta,
                   a_sigma, b_sigma, mu_gamma, sigma_gamma, mu_z, sigma_z, mu_w, sigma_w,
                   lambda0, lambda1, beta, theta, sigma, gamma, z0, z1, w,
                   I, N, G, mlen, mseg, mH, mY);

      // x.format(CommaInitFmt);
      osample << chain_id << ", " << lambda0.format(CSVFormat) << ", " << lambda1.format(CSVFormat)
              << ", " << theta.format(CSVFormat)
              << ", " << beta.format(CSVFormat)
              << ", " << z0.format(CSVFormat) << ", " << z1.format(CSVFormat)
              << ", " << w.format(CSVFormat)
              << ", " << gamma.format(CSVFormat)
              << ", " << sigma
              << ", " << lp_
              << std::endl;
    }
  }
  osample.close();

  // writing (overall) acceptance rate... (last column for sigma's acceptance rate = 1)
  osummary.open(fsummary.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }
  // if (is_empty(osummary)) osummary << ".chain, " << "accept_stat\n";
  osummary << chain_id << ", " << acc_lambda0.mean() / (double)num_warmup
           << ", " << acc_lambda1.mean() / (double)num_warmup
           << ", " << acc_beta.mean() / (double)num_warmup
           << ", " << acc_theta.mean() / (double)num_warmup
           << ", " << acc_z0.mean() / (double)num_warmup
           << ", " << acc_z1.mean() / (double)num_warmup
           << ", " << acc_w.mean() / (double)num_warmup
           << ", " << acc_gamma.mean() / (double)num_warmup
           << ", " << acc_sigma
           << std::endl;
  osummary.close();

  return 1;
}
