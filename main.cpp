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

// see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
const static Eigen::IOFormat CSVFormatN(Eigen::StreamPrecision,
                                        Eigen::DontAlignCols, ", ", "\n");

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, ", ", ", ");

const static Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
                                          Eigen::DontAlignCols, ", ", ", ", "",
                                          "", "", "\n");


int main(int argc, const char *argv[]) {

  int num_samples = atoi(argv[1]);
  int num_warmup = atoi(argv[2]);
  int thin = atoi(argv[3]);
  // double my_eps = atof(argv[4]);
  // int my_L = atoi(argv[5]);
  // double min_E = atof(argv[6]);
  // double max_E = atof(argv[7]);

  int num_iter = num_samples + num_warmup;
  int num_print = num_iter / 100;
  int v = 1; // chain ID

  std::ofstream osummary, osample;
  std::stringstream fstart, fsummary, fsample, fwarmup;

  fstart << "./output/start.csv";
  fwarmup << "./output/warmup.csv";
  fsummary << "./ouput/summary.csv";
  fsample << "./output/sample.csv";

  // use current time as seed for random generator
  int rseed = (unsigned int)time(0) / 2;
  boost::ecuyer1988 rng = create_rng(rseed, v);
  (void)rng; // suppress unused var warning

  std::cout << std::fixed << std::setprecision(1);
  osummary.open(fstart.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }
  // if (is_empty(osummary)) inp << ".chain, " << "seed, " << "iter, " << "warm,
  // " << "thin" << std::endl;
  osummary << v << ", " << rseed << ", " << num_samples << ", " << num_warmup
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
  acc_sigma = 0;
  acc_z0.setZero();
  acc_z1.setZero();
  acc_w.setZero();
  acc_gamma.setZero();
  acc_lambda0.setZero();
  acc_lambda1.setZero();

  // // proprosed parameters
  // Eigen::MatrixXd theta_s(N,2);
  // Eigen::MatrixXd beta_s(I,2);
  // double sigma_s;
  // Eigen::MatrixXd z0_s(N,2);
  // Eigen::MatrixXd z1_s(N,2);
  // Eigen::MatrixXd w_s(I,2);
  // Eigen::VectorXd gamma_s(2);
  // Eigen::MatrixXd lambda0_s(I,G);
  // Eigen::MatrixXd lambda1_s(I,G);

  // // log acceptance ratio
  // Eigen::MatrixXd logr_lambda0(I,G);
  // Eigen::MatrixXd logr_lambda1(I,G);

  osummary.open(fwarmup.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }
  // if (is_empty(osummary)) osummary << ".chain, " << "accept_stat\n";
  osummary << v << ", " << 0 << std::endl;
  osummary.close();

  std::cout << "Starting Sampling...\n";
  for (int ii = 1; ii <= num_warmup; ii++) {

    std::clock_t c_start = std::clock();

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

    for (int k = 0; k < N; k++) {

    update_theta( theta(k,0), acc_theta(k,0), sigma,
                  lambda0, beta.col(0), gamma(0), z0.row(k), w,
                  I, mlen, mseg.col(k), mH.col(k), mY.col(k), 0, rng);

    update_theta( theta(k,1), acc_theta(k,1), sigma,
                  lambda1, beta.col(1), gamma(1), z1.row(k), w,
                  I, mlen, mseg.col(k), mH.col(k), mY.col(k), 1, rng);

    }

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

  std::cout << std::fixed << std::setprecision(1);
  osample.open(fsample.str(), std::ios::app);
  if (!osample.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }

  // init: acceptance ratio
  acc_theta.setZero();
  acc_beta.setZero();
  acc_sigma = 0;
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
      // x.format(CommaInitFmt);
      osample << v << ", " << lambda0.format(CSVFormat) << ", " << lambda1.format(CSVFormat)
              << ", " << theta.format(CSVFormat)
              << ", " << beta.format(CSVFormat)
              << ", " << z0.format(CSVFormat) << ", " << z1.format(CSVFormat)
              << std::endl;
    }
  }
  osample.close();
}
