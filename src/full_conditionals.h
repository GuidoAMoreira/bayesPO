#ifndef __FULL_CONDITIONALS_BAYESPO_H__
#define __FULL_CONDITIONALS_BAYESPO_H__

#include <RcppEigen.h>

// Err by excess precision
#define LOG_SQRT_2 0.3465735902799726547086160607290882840374

// Necessary pre-declarations
class Covariance;

// Abstract classes
class LambdaStar {
public:
  double l;

  virtual double update(long ns, double area) = 0;
}; // LambdaStar

class BetaDelta {
public:
  long s;
  Eigen::VectorXd effects;
  bool initializedGP;

  virtual double update(const Eigen::MatrixXd&,
                        const Eigen::MatrixXd&) = 0;
  virtual Eigen::MatrixXd link(const Eigen::MatrixXd&, bool) = 0;
  virtual double updateGPx(const Eigen::MatrixXd&) = 0;
  virtual double updateGP(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                          const Eigen::MatrixXd&) = 0;

  // Non virtual functions
  BetaDelta(Eigen::MatrixXd, Covariance*);
  double updateCovariance();

  // Gaussian Process members
  const Eigen::MatrixXd xPositions;
  Eigen::MatrixXd xDists;
  Covariance *covariance;
  Eigen::MatrixXd Positions;
  Eigen::VectorXd xMarks, marks, offset;
}; // BetaDelta

class Covariance {
public:
  Covariance(bool);
  double sigma2;

  bool used;
  virtual double update(BetaDelta*) {return 0.;} // Not an abstract class
  virtual double calcCov(double) {return 0.;}    // in case there is no GP.
  Eigen::MatrixXd calcCovMatrix(const Eigen::MatrixXd&, const Eigen::MatrixXd&);
  Eigen::MatrixXd calcCovMatrixPoints(const Eigen::MatrixXd&);
  Eigen::MatrixXd applyCovariance(const Eigen::MatrixXd&);
}; // Covariance

// Derived classes
class gamma_prior : public LambdaStar {
  // Prior parameters
  double shape, rate;
public:
  gamma_prior(Rcpp::List);

  double update(long ns, double area);
}; // gamma_prior

class logit_normal : public BetaDelta {
  // Prior parameters
  Eigen::VectorXd mu, Bb;
  Eigen::MatrixXd Sigma;

  Eigen::VectorXd polyagamma; // Used in sampling
public:
  double update(const Eigen::MatrixXd& onesCov,
                const Eigen::MatrixXd& ZerosCov);
  Eigen::MatrixXd link(const Eigen::MatrixXd&, bool);
  double updateGPx(const Eigen::MatrixXd&);
  double updateGP(const Eigen::MatrixXd &pos, const Eigen::MatrixXd &onesC,
                  const Eigen::MatrixXd &zerosC);

  logit_normal(Rcpp::List, Eigen::MatrixXd, Covariance*);
}; // logit_normal

class Exponential_Covariance : public Covariance {
  // Prior parameters
  double shape, rate;
public:
  double phi;
  double update(BetaDelta*);
  double calcCov(double);
}; // Exponential_Covariance

class Gaussian_Covariance : public Covariance {
  double shape, rate;
public:
  double phi;
  double update(BetaDelta*);
  double calcCov(double);
}; // Gaussian_Covariance

class Matern_Covariance : public Covariance {
public:
  double rho, nu;
  double update(BetaDelta*);
  double calcCov(double);
}; // Matern_Covariance

#endif
