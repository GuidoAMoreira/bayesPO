#include "link_functions.h"

double logistic(vec b, vec cov, bool c)
{
  double out = dot(b,join_vert(ones(1),cov)); // Faster than manual for()
  out = -std::log( 1 + std::exp( pow(-1,!c) * out));

  return out;
}
