#include <stan/math.hpp>
#include <iostream>
#include "my_header.h"
// basic file operations
#include <fstream>
#include <cstdlib>  // to process main(arguemnts)
#include <ctime>    // std::time
#include <typeinfo> // typeid(a).name()

// #include "is_empty.h"
#include "create_rng.hpp"
#include "readCSV.h"
  
// // see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
// const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
//                                        Eigen::DontAlignCols, ", ", "\n");

// const static Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
//                                           Eigen::DontAlignCols, ", ", ", ", "",
//                                           "", "", "\n");

int main(int argc, const char *argv[]) {
  
  int iter = atoi(argv[1]);
  int warm = atoi(argv[2]);
  int thin = atoi(argv[3]);
  // double my_eps = atof(argv[4]);
  // int my_L = atoi(argv[5]);
  // double min_E = atof(argv[6]);
  // double max_E = atof(argv[7]);

  int v = 1;
  
  std::ofstream osummary, osample;
  std::stringstream fstart;
  std::stringstream fsummary;
  std::stringstream fsample;

  fstart << "./output/start.csv";
  fsummary << "./ouput/summary.csv";
  fsample << "./output/sample.csv";

  std::cout << "log normal(1 | 2, 3)=" << stan::math::normal_log(1, 2, 3) << std::endl;

  // use current time as seed for random generator
  int rseed = (unsigned int)time(0) / 2;
  boost::ecuyer1988 rng = create_rng(rseed, 0);
  (void)rng; // suppress unused var warning

  std::cout << "Draw a standard normal random variable.." << stan::math::std_normal_rng(rng) << std::endl;

    // cout.precision(1);
  std::cout << std::fixed << std::setprecision(1);
  osummary.open(fstart.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }

  // if (is_empty(osummary)) inp << ".chain, " << "seed, " << "iter, " << "warm, " << "thin" << std::endl;
  osummary << v << ", " << rseed << ", " << iter << ", " << warm << ", " << thin << std::endl;
  osummary.close();
  
  std::cout << "10 iterations would take " << 4.69 << " seconds.\n" << "Adjust your expectations accordingly!\n";

  std::cout << "Iteration: " << std::setw(7) << std::right << 501 << " / "
            << std::setw(7) << std::right << 20000 << " [" << std::setw(3) << 3 << "%]"
            << "  (Warmup)" << std::endl;
  std::cout << "Iteration: " << std::setw(7) << std::right << 2001 << " / "
            << std::setw(7) << std::right << 20000 << " [" << std::setw(3) << 10 << "%]"
            << "  (Sampling)" << std::endl;
    
  osummary.open(fstart.str(), std::ios::app);
  if (!osummary.is_open()) {
    std::cout << "cannot open the log file\n";
    return 0;
  }

  // if (is_empty(osummary)) osummary << ".chain, " << "accept_stat\n";
  osummary << v << ", " << 0 << std::endl;
  osummary.close();

  std::clock_t c_start = std::clock();
  std::clock_t c_end = std::clock();
  double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;

  Eigen::MatrixXd mvar(1,4);
  mvar = readCSV("mvar.csv", 0, 1, 4);
  int I = mvar(0,0);
  int N = mvar(0,1);
  int C = mvar(0,2);
  int G = mvar(0,3);

  Eigen::VectorXd mlen(G-1);
  Eigen::MatrixXd mseg(I,N);
  Eigen::MatrixXd mh(I,N);
  Eigen::MatrixXd mt(I,N);
  Eigen::MatrixXd mi(I,N);
  
  // mlen = readCSV("mlen.csv", G-1, 1);
  // mseg = readCSV("mseg.csv", I, N);
  // mh = readCSV("mh.csv", I, N);
  // mt = readCSV("mt.csv", I, N);
  // mi = readCSV("mi.csv", I, N);

  // mseg = mseg.cast <int> ();  
  std::cout << I << std::endl;
}
