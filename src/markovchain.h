#ifndef __MARKOVCHAIN_BAYESPO_H__
#define __MARKOVCHAIN_BAYESPO_H__

#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif
#include <RcppArmadillo.h>
#include "covariates.h"
using namespace Rcpp;
using arma::vec;
using arma::mat;

// Generic Markov Chain step class
// Parametric form for q(.) and p(.) are attached before MC procedure begins
class mcStep
{
public:
  //// Unchanging attributes
  const R_xlen_t nb, nd;
  const double area;
  const std::vector<R_xlen_t> X;
  const mat zX;

  //// Attributes
  vec beta, delta;
  std::vector<R_xlen_t> U, Xprime; // vectors of indexes
  double lambdaStar, logPosterior;
  //mat zX, zU, wX, wXp, zXXp;
  mat zXXp, wXp, zU, wX;
  unsigned int iteration;
  retrievCovs *background;

  // Initialization
  vec ibeta, idelta;
  double ilambdaStar;

  // Prior values
  List lambdaPars, betaPars, deltaPars; // Covariance matrixes must be inverted

  // Constructor
  mcStep(NumericVector b, NumericVector d, double l, retrievCovs *bb,
         double a, std::vector<R_xlen_t> x, mat zx, mat wx);

  // Methods
  void update();

  // Attaching varying functions depending on model definition
  void attachLambdaFC(double FCl(double *l, R_xlen_t x, R_xlen_t xp, R_xlen_t u, List p, double a));
  void attachBetaFC(double FCb(vec *b, mat *oc, mat *zc, List p, double link(vec b, vec cov, bool c)));
  void attachDeltaFC(double FCd(vec *d, mat *oc, mat *zc, List p, double link(vec b, vec cov, bool c)));
  void attachCalcq(double cq(vec b, vec co, bool c));
  void attachCalcp(double cp(vec d, vec co, bool c));

protected:
  void applyTransitionKernel(); // Markov chain transition kernel

  // Full conditionals
  double FullConditional_processes();
  // Pointers since they are set by the attach functions
  double (*FullConditional_lambda)(double *lambda, R_xlen_t nx, R_xlen_t nxp,
          R_xlen_t nu, List pars, double areaD);
  double (*FullConditional_beta)(vec *beta, mat *OnesCovs, mat *ZeroesCovs,
          List parameters, double link(vec b, vec cov, bool c));
  double (*FullConditional_delta)(vec *delta, mat *OnesCovs, mat *ZeroesCovs,
          List parameters, double link(vec b, vec cov, bool c));

  // Actual link functions. Used by the latent processes generation
  // Pointers since they are set by the attach functions
  double (*calc_q)(vec beta, vec covs, bool complementary);
  double (*calc_p)(vec beta, vec covs, bool complementary);
};

#endif
