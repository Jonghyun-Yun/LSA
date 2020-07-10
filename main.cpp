#include <stan/math.hpp>
#include "tbb/blocked_range2d.h"
#include "tbb/blocked_range3d.h"
#include <iostream>
// basic file operations
#include <fstream>
// #include <cmath>
#include <cstdlib>  // to process main(arguemnts)
#include <ctime>    // std::time
#include <typeinfo> // typeid(a).name()
#include <string>
#include <vector>

// #include "is_empty.h"
#include "create_rng.hpp"
#include "readCSV.h"
#include "update_lambda.h"
#include "update_theta.h"
#include "update_beta.h"
#include "update_z.h"
#include "update_single_z.h"
#include "update_w.h"
#include "update_gamma.h"
#include "update_single_gamma.h"
#include "fun_lp.h"
#include "par_fun_lp.h"

// typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
// typedef Eigen::Map<Eigen::VectorXd> MapVecd;
// typedef Eigen::Block<Eigen::MatrixXd> BlkMatd;

// see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
const static Eigen::IOFormat CSVFormatN(Eigen::StreamPrecision,
                                        Eigen::DontAlignCols, ", ", "\n");

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, ", ", ", ");

const static Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
                                          Eigen::DontAlignCols, ", ", ", ", "",
                                          "", "", "\n");

inline double update_sigma(const double& a_sigma, const double &b_sigma,
                           const Eigen::MatrixXd &theta, const int &N, boost::ecuyer1988 &rng) {
  return stan::math::sqrt(stan::math::inv_gamma_rng(a_sigma + 0.5 * N, b_sigma + 0.5 * theta.array().square().sum(), rng));
}

int main(int argc, const char *argv[]) {

  std::vector<std::string> sarg;
  sarg.assign(argv, argv + argc);

  // commandline argument:

  // use Intel TBB?
  bool RUN_PAR;
  if (sarg[1] == "parallel") {
    RUN_PAR = true;
  } else if (sarg[1] == "serial") {
    RUN_PAR = false;
  } else {
    std::cout << "invalid arguemnt for RUN_PAR.\n" << std::endl;
    return 0;
  }

  // single z or double z
  bool SINGLE_Z;
  if (sarg[2] == "single_z") {
    SINGLE_Z = true;
  } else if (sarg[2] == "double_z") {
    SINGLE_Z = false;
  } else {
    std::cout << "invalid arguemnt for SINGLE_Z.\n" << std::endl;
    return 0;
  }

  // missing?
  bool FULL_OBS;
  if (sarg[3] == "full") {
    FULL_OBS = true;
  } else if (sarg[3] == "sparse") {
    FULL_OBS = false;
  } else {
    std::cout << "invalid arguemnt for FULL_OBS.\n" << std::endl;
    return 0;
  }

  // latent space or not
  bool UPDATE_LATENT;
  if (sarg[4] == "latent") {
    UPDATE_LATENT = true;
  } else if (sarg[4] == "no_latent") {
    UPDATE_LATENT = false;
  } else {
    std::cout << "invalid arguemnt for UPDATE_LATENT.\n" << std::endl;
    return 0;
  }

  // gamma or not
  bool UPDATE_GAMMA;
  if (sarg[5] == "gamma") {
    UPDATE_GAMMA = true;
  } else if (sarg[5] == "no_gamma") {
    UPDATE_GAMMA = false;
  } else {
    std::cout << "invalid arguemnt for UPDATE_GAMMA.\n" << std::endl;
    return 0;
  }

  int chain_id = atoi(argv[6]);
  int num_samples = atoi(argv[7]);
  int num_warmup = atoi(argv[8]);
  int thin = atoi(argv[9]);
  // double my_eps = atof(argv[4]);
  // int my_L = atoi(argv[5]);
  // double min_E = atof(argv[6]);
  // double max_E = atof(argv[7]);

  double lp_;
  int num_print;

  int num_iter = num_samples + num_warmup;
  if (num_iter < 100) {
    num_print = 1;
  }
  else {
    num_print = num_iter / 100;
  }

  std::ofstream osummary, osample;
  std::stringstream fstart, fsummary, fsample;
  // std::streingtream fwarmup;

  fstart << "./output/start.csv";
  // fwarmup << "./output/warmup.csv";
  fsummary << "./output/summary.csv";
  fsample << "./output/sample_chain" << chain_id << ".csv";

  // use current time as seed for random generator
  int rseed = (unsigned int)time(0) / 2;
  boost::ecuyer1988 rng = create_rng(rseed, chain_id);
  (void)rng; // suppress unused var warning

  std::cout << std::fixed << std::setprecision(1);
  osummary.open(fstart.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file:" << fstart.str() << "\n";
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
  // Eigen::MatrixXd mt(I, N);
  Eigen::MatrixXi mY(I, N);
  Eigen::MatrixXi mNA(I, N);

  mlen = readCSV("mlen.csv", G, 1);
  tmp_mseg = readCSV("mseg.csv", I, N);
  mH = readCSV("mh.csv", I, N);
  // mt = readCSV("mt.csv", I, N);
  tmp_mY = readCSV("mi.csv", I, N);

  mseg = tmp_mseg.cast<int>();
  mY = tmp_mY.cast<int>();

  if (FULL_OBS) {
    mNA.setOnes();
  }
  else {
  Eigen::MatrixXd tmp_mNA(I, N);
  tmp_mNA = readCSV("mNA.csv", I, N);
  mNA = tmp_mNA.cast<int>();
  }

  // std::cout << I << ", " << N << ", " << ", " << C << ", " << G
  //           << std::endl;

  // std::cout << mseg.format(CSVFormat)
  //           << std::endl;

  // std::cout << mH.format(CSVFormat)
  //           << std::endl;

  // std::cout << mY.format(CSVFormat)
  //           << std::endl;

  std::cout << "Reading hyperparameters...\n";
  Eigen::MatrixXd tmp_lambda(3*I, G);
  Eigen::MatrixXd tmp_beta(3*I, 2);
  Eigen::MatrixXd tmp_theta(3*N, 2);
  Eigen::VectorXd tmp_sigma(2);
  Eigen::MatrixXd tmp_gamma(3, 2);
  Eigen::MatrixXd tmp_z(3*N, 2);
  Eigen::MatrixXd tmp_w(3*I, 2);

  // std::cout << "Reading hyperparameters for lambda...\n";
  tmp_lambda = readCSV("pj_lambda.csv", 3*I, G);
  // std::cout << "Reading hyperparameters for beta...\n";
  tmp_beta = readCSV("pj_beta.csv", 3*I, 2);
  // std::cout << "Reading hyperparameters for theta...\n";
  tmp_theta = readCSV("pj_theta.csv", 3*N, 2);
  // std::cout << "Reading hyperparameters for sigma...\n";
  tmp_sigma = readCSV("pj_sigma.csv", 2, 1);
  // std::cout << "Reading hyperparameters for gamma...\n";
  tmp_gamma = readCSV("pj_gamma.csv", 3, 2);
  // std::cout << "Reading hyperparameters for z...\n";
  tmp_z = readCSV("pj_z.csv", 3*N, 2);
  // std::cout << "Reading hyperparameters for w...\n";
  tmp_w = readCSV("pj_w.csv", 3*I, 2);

  // std::cout << "Reading hyperparameters for lambda...\n";
  Eigen::MatrixXd a_lambda(I,G);
  Eigen::MatrixXd b_lambda(I,G);
  Eigen::MatrixXd jump_lambda(I,G);

  a_lambda = tmp_lambda.topRows(I);
  b_lambda = tmp_lambda.block(I,0,I,G);
  jump_lambda = tmp_lambda.bottomRows(I);

  // std::cout << "Reading hyperparameters for beta...\n";
  Eigen::MatrixXd mu_beta(I,2);
  Eigen::MatrixXd sigma_beta(I,2);
  Eigen::MatrixXd jump_beta(I,2);

  mu_beta = tmp_beta.topRows(I);
  sigma_beta = tmp_beta.block(I,0,I,2);
  jump_beta = tmp_beta.bottomRows(I);

  // std::cout << "Reading hyperparameters for theta...\n";
  Eigen::MatrixXd mu_theta(N,2);
  Eigen::MatrixXd sigma_theta(N,2);
  Eigen::MatrixXd jump_theta(N,2);

  mu_theta = tmp_theta.topRows(N);
  sigma_theta = tmp_theta.block(N,0,N,2);
  jump_theta = tmp_theta.bottomRows(N);
  
  // std::cout << "Reading hyperparameters for sigma...\n";
  double a_sigma = tmp_sigma(0);
  double b_sigma = tmp_sigma(1);

  // std::cout << "Reading hyperparameters for gamma...\n";
  Eigen::VectorXd mu_gamma(2);
  Eigen::VectorXd sigma_gamma(2);
  Eigen::VectorXd jump_gamma(2);

  mu_gamma = tmp_gamma.row(0);
  sigma_gamma = tmp_gamma.row(1);
  jump_gamma = tmp_gamma.row(2);

  // std::cout << "Reading hyperparameters for z...\n";
  Eigen::MatrixXd mu_z(N,2);
  Eigen::MatrixXd sigma_z(N,2);
  Eigen::MatrixXd jump_z(N,2);

  mu_z = tmp_z.topRows(N);
  sigma_z = tmp_z.block(N,0,N,2);
  jump_z = tmp_z.bottomRows(N);
  
  // std::cout << "Reading hyperparameters for w...\n";
  Eigen::MatrixXd mu_w(I,2);
  Eigen::MatrixXd sigma_w(I,2);
  Eigen::MatrixXd jump_w(I,2);
  
  mu_w = tmp_w.topRows(I);
  sigma_w = tmp_w.block(I,0,I,2);
  jump_w = tmp_w.bottomRows(I);
  
  // std::cout << "Initializing parameters...\n";
  // parameters
  Eigen::MatrixXd theta(N, 2);
  Eigen::MatrixXd beta(I, 2);
  double sigma;
  // Eigen::MatrixXd z0(N, 2);
  // Eigen::MatrixXd z1(N, 2);
  Eigen::MatrixXd w(I, 2);
  Eigen::VectorXd gamma(2);
  // Eigen::MatrixXd lambda0(I, G);
  // Eigen::MatrixXd lambda1(I, G);
  Eigen::MatrixXd lambda(2 * I, G);
  Eigen::MatrixXd cum_lambda(2 * I, N);

  // Initialization
  theta.setZero();
  beta.setZero();
  sigma = 1.0;
  // z0.setZero();
  // z1.setZero();
  w.setZero();
  // lambda0.setOnes();
  // lambda1.setOnes();
  lambda.setOnes();

  // accept_stat
  Eigen::MatrixXd acc_theta = Eigen::MatrixXd::Zero(N, 2);
  Eigen::MatrixXd acc_beta = Eigen::MatrixXd::Zero(I, 2);
  double acc_sigma = 1.0; //sigma from its full conditional
  // Eigen::VectorXd acc_z0 = Eigen::VectorXd::Zero(N);
  // Eigen::VectorXd acc_z1 = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd acc_w = Eigen::VectorXd::Zero(I);
  Eigen::VectorXd acc_gamma;
  if (SINGLE_Z) {
    acc_gamma = Eigen::VectorXd::Zero(1);
  } else {
    acc_gamma = Eigen::VectorXd::Zero(2);
  }
  // Eigen::MatrixXd acc_lambda0 = Eigen::MatrixXd::Zero(I, G);
  // Eigen::MatrixXd acc_lambda1 = Eigen::MatrixXd::Zero(I, G);
  Eigen::MatrixXd acc_lambda = Eigen::MatrixXd::Zero(2*I, G);

  // Eigen::VectorXd acc_all(4*N + 4*I + 2*I*G + 3);
  
  // acc_theta.setZero();
  // acc_beta.setZero();
  // acc_sigma = 1;
  // acc_z0.setZero();
  // acc_z1.setZero();
  // acc_z.setZero();
  // acc_w.setZero();
  // acc_gamma.setZero();
  // acc_lambda0.setZero();
  // acc_lambda1.setZero();
  // acc_lambda.setZero();

  Eigen::MatrixXd z = Eigen::MatrixXd::Zero(2*N, 2);
  Eigen::VectorXd acc_z;
  if (SINGLE_Z) {
    acc_z = Eigen::VectorXd::Zero(N);
    gamma(0) = -1.0;
    gamma(1) = 1.0;
  }
  else {
    gamma.setOnes();
  }

  std::cout << std::fixed << std::setprecision(1);
  osample.open(fsample.str(), std::ios::app);
  if (!osample.is_open()) {
    std::cout << "cannot open the log file:" << fsample.str() << "\n";
    return 0;
  }

  std::cout << "Starting Sampling...\n";
  for (int ii = 1; ii <= num_iter; ii++) {

    std::clock_t c_start = std::clock();

    if (RUN_PAR) {

      // updating lambda...
      tbb::parallel_for(
        tbb::blocked_range3d<int>(0, 2, 0, I, 0, G),
        [&](tbb::blocked_range3d<int> r)
        {
          for (int c=r.pages().begin(); c<r.pages().end(); ++c)
          {
            for (int i=r.rows().begin(); i<r.rows().end(); ++i)
            {
              for (int g=r.cols().begin(); g<r.cols().end(); ++g)
              {

                update_lambda(lambda((c * I) + i, g), acc_lambda((c * I) + i, g),
                              a_lambda(i,g), b_lambda(i,g), jump_lambda(i,g), g,
                              beta(i, c), theta.col(c), gamma(c), z.block(c*N,0,N,2), w.row(i),
                              N, mNA.row(i), mlen(g), mseg.row(i), mH.row(i), mY.row(i), c, rng);

              }
            }
          }
        });


      // updating beta...
      tbb::parallel_for(
        tbb::blocked_range2d<int>(0,I,0,2),
        [&](tbb::blocked_range2d<int> r)
        {
          for (int i=r.rows().begin(); i<r.rows().end(); ++i)
          {
            for (int c=r.cols().begin(); c<r.cols().end(); ++c)
            {

              cum_lambda.row(c*I + i) =
                update_beta(beta(i,c), acc_beta(i,c), mu_beta(i,c), sigma_beta(i,c), jump_beta(i,c),
                            lambda.row(c*I + i), theta.col(c), gamma(c), z.block(c*N,0,N,2), w.row(i),
                            N, mNA.row(i), mlen, mseg.row(i), mH.row(i), mY.row(i), c, rng);

            }
          }
        });

      // updating theta...
      tbb::parallel_for(
        tbb::blocked_range2d<int>(0,N,0,2),
        [&](tbb::blocked_range2d<int> r)
        {
          for (int k=r.rows().begin(); k<r.rows().end(); ++k)
          {
            for (int c=r.cols().begin(); c<r.cols().end(); ++c)
            {

              update_theta( theta(k,c), acc_theta(k,c), mu_theta(k,c), jump_theta(k,c), sigma,
                            cum_lambda.block(c*I,k,I,1),
                            beta.col(c), gamma(c), z.row(c*N + k), w,
                            I,  mNA.col(k), mlen, mseg.col(k), mH.col(k), mY.col(k), c, rng);

            }
          }
        });

      if (UPDATE_LATENT) {

        if (UPDATE_GAMMA) {
          // updating gamma...

          if (SINGLE_Z) {
            update_single_gamma(gamma, acc_gamma, mu_gamma, sigma_gamma, jump_gamma,
                                cum_lambda, beta, theta, z, w, I, N,
                                G, mNA, mlen, mseg, mH, mY, RUN_PAR, rng);

          }
          else {
            update_gamma(gamma, acc_gamma, mu_gamma, sigma_gamma, jump_gamma,
                         cum_lambda, beta, theta, z, w, I, N,
                         G, mNA, mlen, mseg, mH, mY, rng);
          }

        }
       
        // updating z...
        if (SINGLE_Z) {

          tbb::parallel_for(
          tbb::blocked_range<int>(0, N),
          [&](tbb::blocked_range<int> r) {
              for (int k = r.begin(); k < r.end(); ++k) {
                  z.row(k) = update_single_z(
                    z.row(k), acc_z(k), mu_z(k), sigma_z(k), jump_z(k),
                    cum_lambda.col(k), beta, theta.row(k), gamma,
                    w, I, mNA.col(k), mlen, mseg.col(k), mH.col(k),
                    mY.col(k), rng);
                }
            });
          // copy z0 to z1
          z.block(N,0,N,2) = z.block(0,0,N,2);

        }
        else
        {
          tbb::parallel_for(
            tbb::blocked_range2d<int>(0, N, 0, 2),
            [&](tbb::blocked_range2d<int> r) {
              for (int k = r.rows().begin(); k < r.rows().end(); ++k) {
                for (int c = r.cols().begin(); c < r.cols().end(); ++c) {

                  z.row(c * N + k) = update_z(
                    z.row(c * N + k), acc_z(c * N + k), mu_z(k), sigma_z(k),
                    jump_z(k), cum_lambda.block(c * I, k, I, 1), beta.col(c),
                    theta(k, c), gamma(c), w, I, mNA.col(k), mlen, mseg.col(k), mH.col(k),
                    mY.col(k), c, rng);
                }

              }
            });
        }

        // updating w...
        tbb::parallel_for(
            tbb::blocked_range<int>(0, I), [&](tbb::blocked_range<int> r) {
              for (int i = r.begin(); i < r.end(); ++i) {

                w.row(i) = update_w(
                    w.row(i), acc_w(i), mu_w(i), sigma_w(i), jump_w(i),
                    cum_lambda.row(i), cum_lambda.row(I + i), beta.row(i), theta, gamma,
                    z, N, G, mNA.row(i), mlen,
                    mseg.row(i), mH.row(i), mY.row(i), rng);

              }
            });

      } // end of latent space updating

    } // end of RUN_PAR
    else
    {
      // updating lambda...
      for (int c = 0; c < 2; c++) {
        for (int i = 0; i < I; i++) {
          for (int g = 0; g < G; g++) {

                update_lambda(lambda((c * I) + i, g), acc_lambda((c * I) + i, g),
                              a_lambda(i,g), b_lambda(i,g), jump_lambda(i,g), g,
                              beta(i, c), theta.col(c), gamma(c), z.block(c*N,0,N,2), w.row(i),
                              N, mNA.row(i), mlen(g), mseg.row(i), mH.row(i), mY.row(i), c, rng);

          }
        }
      }

      // updating beta...
      for (int i = 0; i < I; i++) {
        for (int c = 0; c < 2; c++) {

          cum_lambda.row(c*I + i) =
            update_beta(beta(i,c), acc_beta(i,c), mu_beta(i,c), sigma_beta(i,c), jump_beta(i,c),
                        lambda.row(c*I + i), theta.col(c), gamma(c), z.block(c*N,0,N,2), w.row(i),
                        N, mNA.row(i), mlen, mseg.row(i), mH.row(i), mY.row(i), c, rng);

        }
      }

      // updating theta...
      for (int k = 0; k < N; k++) {
        for (int c = 0; c < 2; c++) {

              update_theta( theta(k,c), acc_theta(k,c), mu_theta(k,c), jump_theta(k,c), sigma,
                            cum_lambda.block(c*I,k,I,1),
                            beta.col(c), gamma(c), z.row(c*N + k), w,
                            I,  mNA.col(k), mlen, mseg.col(k), mH.col(k), mY.col(k), c, rng);

        }
      }

      if (UPDATE_LATENT) {

        if (UPDATE_GAMMA) {
          // updating gamma...

          if (SINGLE_Z) {
            update_single_gamma(gamma, acc_gamma, mu_gamma, sigma_gamma, jump_gamma,
                                cum_lambda, beta, theta, z, w, I, N,
                                G, mNA, mlen, mseg, mH, mY, RUN_PAR, rng);

          }
          else {
            update_gamma(gamma, acc_gamma, mu_gamma, sigma_gamma, jump_gamma,
                         cum_lambda, beta, theta, z, w, I, N,
                         G, mNA, mlen, mseg, mH, mY, rng);
          }

        }

      // updating z...
        if (SINGLE_Z) {

          for (int k = 0; k < N; k++) {
            for (int c = 0; c < 2; c++) {
              z.row(k) = update_single_z(
                z.row(k), acc_z(k), mu_z(k), sigma_z(k),
                jump_z(k),
                cum_lambda.col(k), beta, theta.row(k), gamma,
                w, I, mNA.col(k), mlen, mseg.col(k), mH.col(k),
                mY.col(k), rng);
            }
          }
          // copy z0 to z1
          z.block(N,0,N,2) = z.block(0,0,N,2);
        }
        else {

          for (int k = 0; k < N; k++) {
            for (int c = 0; c < 2; c++) {

              z.row(c * N + k) = update_z(
                z.row(c * N + k), acc_z(c * N + k), mu_z(k), sigma_z(k),
                jump_z(k), cum_lambda.block(c * I, k, I, 1), beta.col(c),
                theta(k, c), gamma(c), w, I, mNA.col(k), mlen, mseg.col(k), mH.col(k),
                mY.col(k), c, rng);

            }
          }

        }

      // updating w...
      for (int i = 0; i < I; i++) {

                w.row(i) = update_w(
                    w.row(i), acc_w(i), mu_w(i), sigma_w(i), jump_w(i),
                    cum_lambda.row(i), cum_lambda.row(I + i), beta.row(i), theta, gamma,
                    z, N, G, mNA.row(i), mlen,
                    mseg.row(i), mH.row(i), mY.row(i), rng);

      }

      } // end of latent space updating

    }

    // updating sigma
    sigma = update_sigma(a_sigma, b_sigma, theta, N, rng);

    if (ii == 10) {
      std::clock_t c_end = std::clock();
      double time_elapsed_ms = 10.0 * 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
      std::cout << "\n100 iterations would take " << time_elapsed_ms
                << " milliseconds.\n"
                << "Adjust your expectations accordingly!\n\n";
    }
    if (((ii % num_print == 0) && (ii > 10)) || ii ==  num_samples)
      std::cout << "Chain " << std::setw(3) << chain_id << ": "
                << "Iteration: " << std::setw(7) << std::right << ii << " / "
                << std::setw(7) << std::right << num_iter << " ["
                << std::setw(3) << (int)(ii / (double)num_iter * 100) << "%]"
                << (ii<=num_warmup?"  (Warmup)":"  (Sampling)") << std::endl;

    if ((ii % thin == 0) && (ii > num_warmup)) {
// eval log_prob
      if (RUN_PAR) {
        lp_ = par_fun_lp(a_lambda, b_lambda, mu_beta, sigma_beta, mu_theta, sigma_theta,
                         a_sigma, b_sigma, mu_gamma, sigma_gamma, mu_z, sigma_z, mu_w, sigma_w,
                         lambda, cum_lambda, beta, theta, sigma, gamma, z, w,
                         I, N, G, mNA, mlen, mseg, mH, mY, SINGLE_Z, UPDATE_GAMMA);
      }
      else {
        lp_ = fun_lp(a_lambda, b_lambda, mu_beta, sigma_beta, mu_theta, sigma_theta,
                     a_sigma, b_sigma, mu_gamma, sigma_gamma, mu_z, sigma_z, mu_w, sigma_w,
                     lambda, cum_lambda, beta, theta, sigma, gamma, z, w,
                     I, N, G, mNA, mlen, mseg, mH, mY, SINGLE_Z, UPDATE_GAMMA);
      }
      // x.format(CommaInitFmt);
      osample << chain_id << ", " << ii
              << ", " << lambda.format(CSVFormat)
              << ", " << theta.format(CSVFormat)
              << ", " << beta.format(CSVFormat)
              << ", " << z.format(CSVFormat)
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
      std::cout << "cannot open the log file:" << fsummary.str() << "\n";
      return 0;
    }
// if (is_empty(osummary)) osummary << ".chain, " << "accept_stat\n";
  osummary << chain_id
           << ", " << acc_lambda.mean() / (double)num_iter
           << ", " << acc_beta.mean() / (double)num_iter
           << ", " << acc_theta.mean() / (double)num_iter
           << ", " << acc_z.mean() / (double)num_iter
           << ", " << acc_w.mean() / (double)num_iter
           << ", " << acc_gamma.mean() / (double)num_iter
           << ", " << acc_sigma
           << std::endl;
  osummary.close();

  return 1;
}
