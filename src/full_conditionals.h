#ifndef __FULL_CONDITIONALS_BAYESPO_H__
#define __FULL_CONDITIONALS_BAYESPO_H__

#include <RcppEigen.h>

// Err by excess precision
#define LOG_SQRT_2 0.3465735902799726547086160607290882840374

// Abstract classes
class LambdaStar {
public:
  double l;

  virtual double update(int ns, double area) = 0;
}; // LambdaStar

class BetaDelta {
protected:
  Eigen::MatrixXd xPositions;
public:
  int s;
  Eigen::VectorXd effects;

  virtual double update(const Eigen::MatrixXd&,
                        const Eigen::MatrixXd&) = 0;
  virtual Eigen::MatrixXd link(const Eigen::MatrixXd&, bool) = 0;

  // Non virtual functions
  BetaDelta(Eigen::MatrixXd);
}; // BetaDelta

// Derived classes
class gamma_prior : public LambdaStar {
  // Prior parameters
  double shape, rate;
public:
  gamma_prior(Rcpp::List);

  double update(int ns, double area);
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

  logit_normal(Rcpp::List, Eigen::MatrixXd);
}; // logit_normal

#endif
