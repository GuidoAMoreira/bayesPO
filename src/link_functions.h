#ifndef __LINK_FUNCTIONS_BAYESPO_H__
#define __LINK_FUNCTIONS_BAYESPO_H__

/*
 * WARNING: This code has not been thread safe proofed. Do NOT attempt to
 * run with multiple cores.
 */

#include "PolyaGamma.h"
#include <RcppEigen.h>
extern "C" {
#include "safeR.h"
}
#include "samplers.h"
#include "prior.h"

class LinkFunction {
protected:
  Eigen::VectorXd* tempCovariates;
  double candidateLP;
public:
  // Actual link function as a () operator
  virtual double operator()(double lp, bool log = false,
                          bool complementary = false) = 0;
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& lp,
                                   bool log = false,
                                   bool complementary = false) = 0;

  virtual void acceptCandidate(bool success) = 0; // success is the binomial success
  void setCandidate(double lp, Eigen::VectorXd* c) {
    candidateLP = lp;
    tempCovariates = c;
  }

  // This should generate/process the data augmentation for X
  virtual void startup(const Eigen::VectorXd& Xlp, int maxN) = 0;
  // This should sample the parameters. The prior should have the sampler
  virtual Eigen::VectorXd wrapup(Prior* prior) = 0;

  LinkFunction() {}
  virtual ~LinkFunction() = default;
};

class Logit : public LinkFunction {
  PolyaGamma pg;
  const std::vector<Eigen::MatrixXd> xOuterProds;
  const Eigen::VectorXd xMed;

  // Accumulated to be used in the end
  Eigen::MatrixXd V;
  Eigen::VectorXd med;

  /*
   * This is for the constructor's benefit. It stores the outer products of
   * the rows of the observed matrix X, as they'll be used every iteration.
   * If this proves too memory intensive, maybe it should be repeated every
   * iteration? But then again, if it's too memory intensive, it means that
   * repeating every iteration is going to be very computationally intensive.
   */
  std::vector<Eigen::MatrixXd> createOuterProds(const Eigen::MatrixXd& x) {
    std::vector<Eigen::MatrixXd> output(x.rows());
    for (int i = 0; i < x.rows(); i++) {
      output[i] = x.row(i).transpose() * x.row(i);
    }
    return output;
  }
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
    V += *tempCovariates * tempCovariates->transpose() * pg.draw_like_devroye(candidateLP);
    med += (success ? 1. : -1.) * *tempCovariates * 0.5;
  }

  void startup(const Eigen::VectorXd& startLP, int maxN) {
    V = Eigen::MatrixXd::Zero(xOuterProds[0].rows(), xOuterProds[0].rows());
    med = xMed;
    for (int i = 0; i < startLP.size(); i++) {
      V += xOuterProds[i] * pg.draw_like_devroye(startLP(i));
    }
  }
  Eigen::VectorXd wrapup(Prior* prior) {
    return prior->sampleMP(med, V);
  }

  Logit(const Eigen::MatrixXd& x) : pg(PolyaGamma(1)),
  xOuterProds(createOuterProds(x)), xMed(x.colwise().sum() * 0.5) {}
};

class Probit : public LinkFunction {
  Eigen::VectorXd augmentation;
  int currentIdx;

  Eigen::MatrixXd XtX;
  Eigen::VectorXd Xty;

  const Eigen::MatrixXd X;
  const Eigen::MatrixXd XtXog;
public:
  double operator()(double lp, bool log = false, bool complementary = false) {
    return safe_pnorm(lp * (complementary ? -1. : 1.), 0., 1., true, log);
  }
  Eigen::VectorXd operator()(const Eigen::VectorXd& lp, bool log = false,
                           bool complementary = false) {
    Eigen::VectorXd output(lp.size());
    for (int i = 0; i < lp.size(); i++) {
      output(i) = (*this)(lp(i), log, complementary);
    }
    return output;
  }

  void acceptCandidate(bool success) {
    augmentation(currentIdx) = rtnorm(candidateLP, 1., success);
    XtX += *tempCovariates * tempCovariates->transpose();
    Xty += *tempCovariates * augmentation(currentIdx++);
  }

  void startup(const Eigen::VectorXd& startLP, int maxN) {
    augmentation.resize(maxN);
    Xty = Eigen::VectorXd::Zero(XtX.rows(), 1);
    for (currentIdx = 0; currentIdx < startLP.size(); currentIdx++) {
      augmentation(currentIdx) = rtnorm(startLP(currentIdx), 1., true);
      Xty += X.row(currentIdx) * augmentation(currentIdx);
    }
    XtX = XtXog;
  }
  Eigen::VectorXd wrapup(Prior* prior) {
    augmentation.resize(currentIdx);
    return prior->sampleMP(XtX.llt().solve(Xty), XtX);
  }

  Probit(const Eigen::MatrixXd& x) : LinkFunction(), X(x),
  XtXog(x.transpose() * x) {}
};

#endif
