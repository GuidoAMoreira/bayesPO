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
  const std::vector<int> selInt, selObs;
  int ncell, nvar;

  // Virtual methods
  virtual Eigen::VectorXd retrieveInt(int) = 0;
  virtual Eigen::VectorXd retrieveObs(int) = 0;
  void putInt(Eigen::MatrixXd&, std::vector<int>&, int, int, int);
  void putObs(Eigen::MatrixXd&, std::vector<int>&, int, int, int);

  // Universal methods
  virtual int pickRandomPoint();
  virtual Eigen::VectorXi pickRandomPoint(int n);
  void addAcceptedXprime(int);
  void addAcceptedXprime(Eigen::VectorXi points) {
    for (int i = 0; i < points.size(); i++)
      addAcceptedXprime(points(i));
  }
  Eigen::MatrixXd retrieveInt(const Eigen::VectorXi&);
  Eigen::MatrixXd retrieveObs(const Eigen::VectorXi&);
  void resetUnobservedCounts() {
    unObservedCounts.setZero();
  };

  // Constructor
  retrievCovs(std::vector<int> si, std::vector<int> so);
  retrievCovs();

  // getters
  Eigen::VectorXd getUnobservedCounts() {return unObservedCounts;}

protected:
  SEXP covs;
  double *c;
  int nInt;
  int nObs;
  Eigen::VectorXd unObservedCounts;
}; // retrievCovs

// When the covariates are in an integer matrix
class retrievCovs_intMatrix : public retrievCovs
{
public:
  // Methods
  Eigen::VectorXd retrieveInt(int ind);
  Eigen::VectorXd retrieveObs(int ind);

  // Constructor
  retrievCovs_intMatrix(SEXP inp, std::vector<int> si,
                        std::vector<int> so);

private:
  int *c;
}; // retrievCovs_intMatrix

// When the covariates are in an double matrix
class retrievCovs_doubleMatrix : public retrievCovs
{
public:

  // Methods
  Eigen::VectorXd retrieveInt(int ind);
  Eigen::VectorXd retrieveObs(int ind);

  // Constructor
  retrievCovs_doubleMatrix(SEXP inp, std::vector<int> si,
                           std::vector<int> so);
}; // retrievCovs_doubleMatrix

// When the covariates are standard normal
class retrievCovs_normal : public retrievCovs
{
public:
  // Methods
  Eigen::VectorXd retrieveInt(int ind);
  Eigen::VectorXd retrieveObs(int ind);

  // Constructor
  retrievCovs_normal(std::vector<int> si, std::vector<int> so);
}; // retrievCovs_normal

#endif
