#include "covariates.h"

//// Constructors ////
retrievCovs::retrievCovs(std::vector<long> si,
                         std::vector<long> so) : selInt(si),selObs(so) {}
retrievCovs::retrievCovs() {}

retrievCovs_intMatrix::retrievCovs_intMatrix(SEXP inp, std::vector<long> si,
                                             std::vector<long> so) :
  retrievCovs(si,so)
{covs = inp; c = INTEGER(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0]; nvar = INTEGER(dim)[1];}

retrievCovs_doubleMatrix::retrievCovs_doubleMatrix(SEXP inp,
                                                   std::vector<long> si,
                                                   std::vector<long> so) :
  retrievCovs(si,so)
{covs = inp; c = REAL(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0]; nvar = INTEGER(dim)[1];}

retrievCovs_normal::retrievCovs_normal(std::vector<long> si,
                                       std::vector<long> so, long ni,
                                       long no) : retrievCovs(si,so),
                                       n_var_intens(ni),n_var_obs(no) {}

//// Methods ////
//// Base class ////
long retrievCovs::pickRandomPoint()
{return long(R::runif(0,1) * ncell);}

Eigen::VectorXi retrievCovs::pickRandomPoint(long n)
{
  Eigen::VectorXi out(n);
  for (R_xlen_t i = 0; i < n; i++)
    out[i] = pickRandomPoint();

  return out;
}

// Integer Matrix
Eigen::VectorXd retrievCovs_intMatrix::retrieveInt(long ind)
{
  Eigen::VectorXd output(selInt.size());
  for (R_xlen_t i = 0; i < selInt.size(); ++i)
    output[i] = double(c[selInt[i]*ncell + ind]);

  return output;
}

Eigen::VectorXd retrievCovs_intMatrix::retrieveObs(long ind)
{
  Eigen::VectorXd output(selObs.size());
  for (R_xlen_t i = 0; i < selObs.size(); ++i)
    output[i] = double(c[selObs[i]*ncell + ind]);

  return output;
}

// Double Matrix
Eigen::VectorXd retrievCovs_doubleMatrix::retrieveInt(long ind)
{
  Eigen::VectorXd output(selInt.size());
  for (R_xlen_t i = 0; i < selInt.size(); ++i)
    output[i] = c[selInt[i]*ncell + ind];

  return output;
}

Eigen::VectorXd retrievCovs_doubleMatrix::retrieveObs(long ind)
{
  Eigen::VectorXd output(selObs.size());
  for (R_xlen_t i = 0; i < selObs.size(); ++i)
    output[i] = c[selObs[i]*ncell + ind];

  return output;
}

// Normal
Eigen::VectorXd retrievCovs_normal::retrieveInt(long ind)
{
  Eigen::VectorXd output(n_var_intens);
  for (R_xlen_t i = 0; i < n_var_intens; i++)
   output[i] = R::rnorm(0,1);

  return output;
}

Eigen::VectorXd retrievCovs_normal::retrieveObs(long ind)
{
  Eigen::VectorXd output(n_var_obs);
  for (R_xlen_t i = 0; i < n_var_obs; i++)
    output[i] = R::rnorm(0,1);

  return output;
}
