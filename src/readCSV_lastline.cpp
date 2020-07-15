// Copyright Robert Eisele 2017
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <Eigen/Dense>
#include "readCSV_lastline.h"

std::istream& ignoreline(std::ifstream& in, std::ifstream::pos_type& pos)
{
  pos = in.tellg();
  return in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

std::string getLastLine(std::ifstream& in)
{
  std::ifstream::pos_type pos = in.tellg();

  std::ifstream::pos_type lastPos;
  while (in >> std::ws && ignoreline(in, lastPos))
    pos = lastPos;

  in.clear();
  in.seekg(pos);

  std::string line;
  std::getline(in, line);
  return line;
}

Eigen::VectorXd readCSV_lastline(std::string file, int cols) {

  std::ifstream in(file);
  
  std::string line;

  int col = 0;

  Eigen::VectorXd res = Eigen::VectorXd(cols);

  if (in.is_open()) {

    line = getLastLine(in);
    char *ptr = (char *) line.c_str();
    int len = line.length();

    char *start = ptr;
    for (int i = 0; i < len; i++) {

      if (ptr[i] == ',') {
        res(col++) = atof(start);
        start = ptr + i + 1;
      }
    }
    res(col) = atof(start);

    in.close();
  }
  return res;
}

