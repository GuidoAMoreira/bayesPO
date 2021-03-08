#include "covariates.h"
#include "markovchain.h"
#include "full_conditionals.h"
#include "link_functions.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using arma::mat;
using arma::vec;

// Create the mcStep object to be updated. Requires the pointer to a retrievCovs
// variable
mcStep startup_mcmc(NumericVector beta, NumericVector delta, double lambda,
                    std::string b_updater, std::string d_updater,
                    std::string l_updater, List parB, List parD, List parL,
                    std::string xClass, SEXP xValues, IntegerVector xI,
                    IntegerVector xO, retrievCovs *covars, double aD);

// Formating x functions
List importX_double(NumericMatrix x, R_xlen_t nb, R_xlen_t nd, IntegerVector xI, IntegerVector xO);
List importX_int(IntegerMatrix x, R_xlen_t nb, R_xlen_t nd, IntegerVector xI, IntegerVector xO);
List determineX(NumericVector x, R_xlen_t nb, R_xlen_t nd, retrievCovs *b);
