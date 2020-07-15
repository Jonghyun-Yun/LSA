#ifndef __READCSV_LASTLINE_H_
#define __READCSV_LASTLINE_H_

std::istream& ignoreline(std::ifstream& in, std::ifstream::pos_type& pos);

std::string getLastLine(std::ifstream& in);

Eigen::VectorXd readCSV_lastline(std::string file, int cols);

#endif // __READCSV_LASTLINE_H_
