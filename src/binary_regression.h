#ifndef __BINARY_REGRESSION_BAYESPO_H__
#define __BINARY_REGRESSION_BAYESPO_H__

#include <RcppEigen.h>
extern "C" {
#include "safeR.h"
}
#include "utils.h"
#include "PolyaGamma.h"
#include "samplers.h"

class Prior;

class LinkFunction {
protected:
  double candidateLP;
  Eigen::VectorXd augmentation;
  int currentIdx;
public:
  // Actual link function
  virtual double operator()(double lp, bool log = false,
                          bool complementary = false) = 0;
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& lp,
                                   bool log = false,
                                   bool complementary = false) = 0;

  virtual void acceptCandidate(bool success) = 0;
  void setCandidate(double lp) {candidateLP = lp;}

  virtual void startup(const Eigen::VectorXd& X) = 0;

  // Getters
  Eigen::VectorXd getAugmentation() {return augmentation;}
  double getAugmentation(int i) {return augmentation(i);}

  LinkFunction() {}
};

class Logit : public LinkFunction {
  PolyaGamma pg;
  Eigen::MatrixXd V;
public:
  double operator()(double lp, bool log = false, bool complementary = false) {
    double logVal = -log1p(exp(lp * (complementary ? 1. : -1.)));
    if (log) return logVal; else return exp(logVal);
  }
  Eigen::VectorXd operator()(const Eigen::VectorXd& lp, bool log = false,
                           bool complementary = false) {
    Eigen::VectorXd logVals =
      -(lp * (complementary ? 1. : -1.)).array().exp().log1p();
    if (log) return logVals; else return logVals.array().exp();
  }

  void acceptCandidate(bool success) {
    augmentation(currentIdx++) = pg.draw_like_devroye(candidateLP);
  }

  void startup(const Eigen::VectorXd& startLP, int maxN) {
    augmentation.resize(maxN);
    for (currentIdx = 0; currentIdx < startLP.size(); currentIdx++) {
      augmentation(currentIdx) = pg.draw_like_devroye(startLP(currentIdx));
    }
  }

  Logit() : pg(PolyaGamma(1)) {}
};

class Probit : public LinkFunction {
  double operator()(double lp, bool log = false, bool complementary = false) {
    return safe_pnorm(lp * (complementary ? -1. : 1.), 0., 1., true, log);
  }
  Eigen::VectorXd operator()(const Eigen::VectorXd& lp, bool log = false,
                           bool complementary = false) {
    return safe_pnorm(lp * (complementary ? -1. : 1.), 0., 1., true, log);
  }

  void acceptCandidate(bool success) {
    augmentation(currentIdx++) = rtnorm(candidateLP, 1., success);
  }

  void startup(const Eigen::VectorXd& startLP, int maxN) {
    augmentation.resize(maxN);
    for (currentIdx = 0; currentIdx < startLP.size(); currentIdx++) {
      augmentation(currentIdx) = rtnorm(startLP(currentIdx), 1., true);
    }
  }

  Probit() {}
};

class BinaryRegression {
  virtual Eigen::VectorXd getLP(const Eigen::MatrixXd&) = 0; // LP = Linear Predictor
protected:
  Eigen::VectorXd parameters;
  double candidateLP; // LP = Linear Predictor

  Prior* prior;
  LinkFunction* linkFunction;
public:
  /*
   * This function should receive a vector of covariates and return the
   * success probability under the model.
   */
  double considerCandidate(const Eigen::VectorXd& c) {
    double lp = getLP(c)(0);
    linkFunction->setCandidate(lp);
    return (*linkFunction)(lp, true);
  }
  virtual void acceptCandidate(bool success) = 0; // success is the binomial success

  // Constructor
  BinaryRegression() {}
  BinaryRegression(Prior* p, LinkFunction* lf, Eigen::VectorXd init) :
    prior(p), linkFunction(lf), parameters(init) {}
};



#endif
