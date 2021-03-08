#ifndef __RETRIEVE_COVARIATES_BAYESPO_H__
#define __RETRIEVE_COVARIATES_BAYESPO_H__
#include <RcppArmadillo.h>
using arma::vec;
using arma::mat;

/*
To retrieve the covariates regardless of its class, define a generic class with a virtual
retrieve function, then subclasses with the appropriate methods. Declare pointers to the base
class and pass on pointers to the appropriate derived class.
*/
class retrievCovs
{
public:
  // Attributes
  const std::vector<R_xlen_t> selInt, selObs;
  R_xlen_t ncell, nvar;

  // Virtual methods
  virtual vec retrieveInt(R_xlen_t ind) {vec out; return out;} // Unused definition
  virtual vec retrieveObs(R_xlen_t ind) {vec out; return out;} // Unused definition

  // Constructor
  retrievCovs(std::vector<R_xlen_t> si, std::vector<R_xlen_t> so);
  retrievCovs();

protected:
  SEXP covs;
  double *c;
};

// When the covariates are in an integer matrix
class retrievCovs_intMatrix : public retrievCovs
{
public:

  vec retrieveInt(R_xlen_t ind);
  vec retrieveObs(R_xlen_t ind);

  // Constructor
  retrievCovs_intMatrix(SEXP inp, std::vector<R_xlen_t> si, std::vector<R_xlen_t> so);

private:
  int *c;
};

// When the covariates are in an double matrix
class retrievCovs_doubleMatrix : public retrievCovs
{
public:

  vec retrieveInt(R_xlen_t ind);
  vec retrieveObs(R_xlen_t ind);

  // Constructor
  retrievCovs_doubleMatrix(SEXP inp, std::vector<R_xlen_t> si, std::vector<R_xlen_t> so);
};

// When the covariates are standard normal
class retrievCovs_normal : public retrievCovs
{
public:
  // Only sizes matter
  R_xlen_t n_var_intens, n_var_obs;

  vec retrieveInt(R_xlen_t ind);
  vec retrieveObs(R_xlen_t ind);

  // Constructor
  retrievCovs_normal(std::vector<R_xlen_t> si, std::vector<R_xlen_t> so,R_xlen_t ni, R_xlen_t no);
};

#endif
