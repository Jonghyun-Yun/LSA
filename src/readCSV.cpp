// Copyright Robert Eisele 2017
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "readCSV.h"

auto readCSV(std::string file, typename T, std::string type, int rows, int cols) {

  std::ifstream in(file);
  
  std::string line;

  int row = 0;
  int col = 0;

  Eigen::Matrix<T, rows, cols> res;

  if (in.is_open()) {

    while (std::getline(in, line)) {

      char *ptr = (char *) line.c_str();
      int len = line.length();

      col = 0;

      char *start = ptr;
      for (int i = 0; i < len; i++) {

        if (ptr[i] == ',') {
          res(row, col++) = atof(start);
          start = ptr + i + 1;
        }
      }
      res(row, col) = atof(start);

      row++;
    }

    in.close();
  }
  return res;
}

