#include "covariates.h"

//// Constructors ////
retrievCovs::retrievCovs(std::vector<R_xlen_t> si,
                         std::vector<R_xlen_t> so) : selInt(si),selObs(so) {}
retrievCovs::retrievCovs() {}

retrievCovs_intMatrix::retrievCovs_intMatrix(SEXP inp, std::vector<R_xlen_t> si,
                                             std::vector<R_xlen_t> so) : retrievCovs(si,so)
{covs = inp; c = INTEGER(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0];nvar = INTEGER(dim)[1];}

retrievCovs_doubleMatrix::retrievCovs_doubleMatrix(SEXP inp,
                                                   std::vector<R_xlen_t> si,
                                                   std::vector<R_xlen_t> so) : retrievCovs(si,so)
{covs = inp; c = REAL(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0];nvar = INTEGER(dim)[1];}

retrievCovs_normal::retrievCovs_normal(std::vector<R_xlen_t> si,
                                       std::vector<R_xlen_t> so,R_xlen_t ni,
                                       R_xlen_t no) : retrievCovs(si,so),n_var_intens(ni),n_var_obs(no) {}

//// Methods ////
// Integer Matrix
vec retrievCovs_intMatrix::retrieveInt(R_xlen_t ind)
{
  vec output(selInt.size());
  for (long unsigned i = 0;i<selInt.size();++i)
    output[i] = double(c[selInt[i]*ncell + ind]);

  return output;
}

vec retrievCovs_intMatrix::retrieveObs(R_xlen_t ind)
{
  vec output(selObs.size());
  for (long unsigned i = 0;i<selObs.size();++i)
    output[i] = double(c[selObs[i]*ncell + ind]);

  return output;
}

// Double Matrix
vec retrievCovs_doubleMatrix::retrieveInt(R_xlen_t ind)
{
  vec output(selInt.size());
  for (long unsigned i = 0;i<selInt.size();++i)
    output[i] = c[selInt[i]*ncell + ind];

  return output;
}

vec retrievCovs_doubleMatrix::retrieveObs(R_xlen_t ind)
{
  vec output(selObs.size());
  for (long unsigned i = 0;i<selObs.size();++i)
    output[i] = c[selObs[i]*ncell + ind];

  return output;
}

// Normal
vec retrievCovs_normal::retrieveInt(R_xlen_t ind)
{
  vec output(n_var_intens);
  output = arma::randn(n_var_intens);

  return output;
}

vec retrievCovs_normal::retrieveObs(R_xlen_t ind)
{
  vec output(n_var_obs);
  output = arma::randn(n_var_obs);

  return output;
}
