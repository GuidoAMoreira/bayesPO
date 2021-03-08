#ifndef __LINKFUNC_BAYESPO_H__
#define __LINKFUNC_BAYESPO_H__

#include <RcppArmadillo.h>
using namespace arma;

double logistic(vec b, vec cov, bool c); // c is for complementary probability

#endif
