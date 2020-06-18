#include <stan/math.hpp>
#include <iostream>
#include "create_rng.hpp"

int main() {
  std::cout << "log normal(1 | 2, 3)=" << stan::math::normal_log(1, 2, 3) << std::endl;

  // use current time as seed for random generator
  int rseed = (unsigned int)time(0) / 2;

  boost::ecuyer1988 rng = my::create_rng(rseed, 0);
  (void)rng; // suppress unused var warning

  std::cout << "Draw a standard normal random variable.." << stan::math::std_normal_rng(rng) << std::endl;
}
