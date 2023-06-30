#include "covariates.h"

//// Constructors ////
retrievCovs::retrievCovs(std::vector<long> si,
                         std::vector<long> so) : selInt(si),selObs(so),
                         nInt(si.size()), nObs(so.size()){}
retrievCovs::retrievCovs() {}

retrievCovs_intMatrix::retrievCovs_intMatrix(SEXP inp, std::vector<long> si,
                                             std::vector<long> so) :
  retrievCovs(si,so)
{covs = inp; c = INTEGER(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0]; nvar = INTEGER(dim)[1];
 unObservedCounts = Eigen::MatrixXd::Constant(ncell, 1, 0);}

retrievCovs_doubleMatrix::retrievCovs_doubleMatrix(SEXP inp,
                                                   std::vector<long> si,
                                                   std::vector<long> so) :
  retrievCovs(si,so)
{covs = inp; c = REAL(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0]; nvar = INTEGER(dim)[1];
 unObservedCounts = Eigen::MatrixXd::Constant(ncell, 1, 0);}

retrievCovs_normal::retrievCovs_normal(std::vector<long> si,
                                       std::vector<long> so) :
  retrievCovs(si,so) {
    unObservedCounts = Eigen::MatrixXd::Constant(ncell, 1, 0);
  }

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

void retrievCovs::addAcceptedXprime(long point) {
  unObservedCounts(point) += 1;
}

void retrievCovs::putInt(Eigen::MatrixXd& covs, std::vector<long>& ind,
                                         long start, long finish, long startMat)
{
  long starter = startMat - start;
  for (R_xlen_t i = start; i <= finish ; i++)
    covs.row(starter + i) = retrieveInt(ind[i]);
}

void retrievCovs::putObs(Eigen::MatrixXd& covs, std::vector<long>& ind,
                                         long start, long finish, long startMat)
{
  long starter = startMat - start;
  for (R_xlen_t i = start; i <= finish ; i++)
    covs.row(starter + i) = retrieveObs(ind[i]);
}

// Integer Matrix
Eigen::VectorXd retrievCovs_intMatrix::retrieveInt(long ind)
{
  Eigen::VectorXd output(nInt);
  for (R_xlen_t i = 0; i < selInt.size(); ++i)
    output[i] = double(c[selInt[i]*ncell + ind]);

  return output;
}

Eigen::VectorXd retrievCovs_intMatrix::retrieveObs(long ind)
{
  Eigen::VectorXd output(nObs);
  for (R_xlen_t i = 0; i < selObs.size(); ++i)
    output[i] = double(c[selObs[i]*ncell + ind]);

  return output;
}

// Double Matrix
Eigen::VectorXd retrievCovs_doubleMatrix::retrieveInt(long ind)
{
  Eigen::VectorXd output(nInt);
  for (R_xlen_t i = 0; i < selInt.size(); ++i)
    output[i] = c[selInt[i]*ncell + ind];

  return output;
}

Eigen::VectorXd retrievCovs_doubleMatrix::retrieveObs(long ind)
{
  Eigen::VectorXd output(nObs);
  for (R_xlen_t i = 0; i < selObs.size(); ++i)
    output[i] = c[selObs[i]*ncell + ind];

  return output;
}

// Normal
Eigen::VectorXd retrievCovs_normal::retrieveInt(long ind)
{
  Eigen::VectorXd output(nInt);
  for (R_xlen_t i = 0; i < nInt; i++)
   output[i] = R::rnorm(0,1);

  return output;
}

Eigen::VectorXd retrievCovs_normal::retrieveObs(long ind)
{
  Eigen::VectorXd output(nObs);
  for (R_xlen_t i = 0; i < nObs; i++)
    output[i] = R::rnorm(0,1);

  return output;
}
