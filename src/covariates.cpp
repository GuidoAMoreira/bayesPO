#include "covariates.h"
#include "samplers.h"

//// Constructors ////
retrievCovs::retrievCovs(std::vector<int> si,
                         std::vector<int> so) : selInt(si),selObs(so),
                         nInt(si.size()), nObs(so.size()){}
retrievCovs::retrievCovs() {}

retrievCovs_intMatrix::retrievCovs_intMatrix(SEXP inp, std::vector<int> si,
                                             std::vector<int> so) :
  retrievCovs(si,so)
{covs = inp; c = INTEGER(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0]; nvar = INTEGER(dim)[1];
 unObservedCounts = Eigen::MatrixXd::Zero(ncell, 1);}

retrievCovs_doubleMatrix::retrievCovs_doubleMatrix(SEXP inp,
                                                   std::vector<int> si,
                                                   std::vector<int> so) :
  retrievCovs(si,so)
{covs = inp; c = REAL(covs); SEXP dim = Rf_getAttrib( inp, R_DimSymbol ) ;
 ncell = INTEGER(dim)[0]; nvar = INTEGER(dim)[1];
 unObservedCounts = Eigen::MatrixXd::Zero(ncell, 1);}

retrievCovs_normal::retrievCovs_normal(std::vector<int> si,
                                       std::vector<int> so) :
  retrievCovs(si,so) {
    unObservedCounts = Eigen::MatrixXd::Zero(ncell, 1);
  }

Eigen::MatrixXd retrievCovs::retrieveInt(const Eigen::VectorXi& indices) {
  Eigen::MatrixXd out(indices.size(), selInt.size());
  for (int i = 0; i < indices.size(); i++)
    out.row(i) = retrieveInt(indices(i)).transpose();

  return out;
}

Eigen::MatrixXd retrievCovs::retrieveObs(const Eigen::VectorXi& indices) {
  Eigen::MatrixXd out(indices.size(), selObs.size());
  for (int i = 0; i < indices.size(); i++)
    out.row(i) = retrieveObs(indices(i)).transpose();

  return out;
}

//// Methods ////
//// Base class ////
int retrievCovs::pickRandomPoint()
{return int(runif(0., 1.) * ncell);}

Eigen::VectorXi retrievCovs::pickRandomPoint(int n)
{
  Eigen::VectorXi out(n);
  for (int i = 0; i < n; i++)
    out(i) = pickRandomPoint();

  return out;
}

void retrievCovs::addAcceptedXprime(int point) {
  unObservedCounts(point) += 1;
}

void retrievCovs::putInt(Eigen::MatrixXd& covs, std::vector<int>& ind,
                                         int start, int finish, int startMat)
{
  int starter = startMat - start;
  for (R_xlen_t i = start; i <= finish ; i++)
    covs.row(starter + i) = retrieveInt(ind[i]);
}

void retrievCovs::putObs(Eigen::MatrixXd& covs, std::vector<int>& ind,
                                         int start, int finish, int startMat)
{
  int starter = startMat - start;
  for (R_xlen_t i = start; i <= finish ; i++)
    covs.row(starter + i) = retrieveObs(ind[i]);
}

// Integer Matrix
Eigen::VectorXd retrievCovs_intMatrix::retrieveInt(int ind)
{
  Eigen::VectorXd output(nInt);
  for (R_xlen_t i = 0; i < selInt.size(); ++i)
    output[i] = double(c[selInt[i]*ncell + ind]);

  return output;
}

Eigen::VectorXd retrievCovs_intMatrix::retrieveObs(int ind)
{
  Eigen::VectorXd output(nObs);
  for (R_xlen_t i = 0; i < selObs.size(); ++i)
    output[i] = double(c[selObs[i]*ncell + ind]);

  return output;
}

// Double Matrix
Eigen::VectorXd retrievCovs_doubleMatrix::retrieveInt(int ind)
{
  Eigen::VectorXd output(nInt);
  for (R_xlen_t i = 0; i < selInt.size(); ++i)
    output[i] = c[selInt[i]*ncell + ind];

  return output;
}

Eigen::VectorXd retrievCovs_doubleMatrix::retrieveObs(int ind)
{
  Eigen::VectorXd output(nObs);
  for (R_xlen_t i = 0; i < selObs.size(); ++i)
    output[i] = c[selObs[i] * ncell + ind];

  return output;
}

// Normal
Eigen::VectorXd retrievCovs_normal::retrieveInt(int ind)
{
  Eigen::VectorXd output(nInt);
  for (R_xlen_t i = 0; i < nInt; i++)
   output[i] = R::rnorm(0,1);

  return output;
}

Eigen::VectorXd retrievCovs_normal::retrieveObs(int ind)
{
  Eigen::VectorXd output(nObs);
  for (R_xlen_t i = 0; i < nObs; i++)
    output[i] = R::rnorm(0,1);

  return output;
}
