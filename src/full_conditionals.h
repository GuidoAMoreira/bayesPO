#ifndef __FULL_CONDITIONALS_BAYESPO_H__
#define __FULL_CONDITIONALS_BAYESPO_H__

#include <RcppArmadillo.h>
#include "PolyaGamma.h"
#include "covariates.h"
using namespace arma;

// Full conditional for Lambda*: l_prior
double l_gamma(double *l, R_xlen_t x, R_xlen_t xp, R_xlen_t u, List p, double aD);

// Full conditional for beta or delta: bd_linkFunction_prior
double bd_logit_normal(vec *bd, mat *oc, mat *zc, List p, double link(vec b, vec cov, bool c));

#endif
