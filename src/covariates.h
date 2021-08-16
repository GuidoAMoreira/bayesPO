#ifndef __RETRIEVE_COVARIATES_BAYESPO_H__
#define __RETRIEVE_COVARIATES_BAYESPO_H__
#include <RcppEigen.h>

/*
To retrieve the covariates regardless of its class, define a generic class with a virtual
retrieve function, then subclasses with the appropriate methods. Declare pointers to the base
class and pass on pointers to the appropriate derived class.
*/
class retrievCovs
{
public:
  // Attributes
  const std::vector<long> selInt, selObs;
  long ncell, nvar;

  // Virtual methods
  virtual Eigen::VectorXd retrieveInt(long) = 0;
  virtual Eigen::VectorXd retrieveObs(long) = 0;
  void putInt(Eigen::MatrixXd&, std::vector<long>&, long, long, long);
  void putObs(Eigen::MatrixXd&, std::vector<long>&, long, long, long);

  // Universal methods
  virtual long pickRandomPoint();
  virtual Eigen::VectorXi pickRandomPoint(long n);

  // Constructor
  retrievCovs(std::vector<long> si, std::vector<long> so);
  retrievCovs();

protected:
  SEXP covs;
  double *c;
  long nInt;
  long nObs;
}; // retrievCovs

// When the covariates are in an integer matrix
class retrievCovs_intMatrix : public retrievCovs
{
public:
  // Methods
  Eigen::VectorXd retrieveInt(long ind);
  Eigen::VectorXd retrieveObs(long ind);

  // Constructor
  retrievCovs_intMatrix(SEXP inp, std::vector<long> si,
                        std::vector<long> so);

private:
  int *c;
}; // retrievCovs_intMatrix

// When the covariates are in an double matrix
class retrievCovs_doubleMatrix : public retrievCovs
{
public:

  // Methods
  Eigen::VectorXd retrieveInt(long ind);
  Eigen::VectorXd retrieveObs(long ind);

  // Constructor
  retrievCovs_doubleMatrix(SEXP inp, std::vector<long> si,
                           std::vector<long> so);
}; // retrievCovs_doubleMatrix

// When the covariates are standard normal
class retrievCovs_normal : public retrievCovs
{
public:
  // Methods
  Eigen::VectorXd retrieveInt(long ind);
  Eigen::VectorXd retrieveObs(long ind);

  // Constructor
  retrievCovs_normal(std::vector<long> si, std::vector<long> so);
}; // retrievCovs_normal

#endif
