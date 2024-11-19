#ifndef __BINARY_REGRESSION_BAYESPO_H__
#define __BINARY_REGRESSION_BAYESPO_H__

/*
 * WARNING: This code has not been thread safe proofed. Do NOT attempt to
 * run with multiple cores.
 */

#include <RcppEigen.h>
#include "utils.h"
#include "link_functions.h"
#include "prior.h"

// This class is abstract because it can be either a (generalized) linear
// regression or, say, a BART.
class BinaryRegression {
  double logPost;
  double transformedLP;
  virtual Eigen::VectorXd getLP(const Eigen::MatrixXd&) = 0; // LP = Linear Predictor

  Prior* prior;
  LinkFunction* linkFunction;
protected:
  Eigen::VectorXd parameters;
  double candidateLP; // LP = Linear Predictor. Not sure why this is not private. Maybe someday I'll know.
public:
  /*
   * These functions deal with receiving covariates, considering (and storing)
   * candidate points, accepting the candidates, then updating the relevant
   * parameters after all candidates have been considered.
   *
   * Candidates are handled differently depending on link function, so it is
   * passed on to them.
   */
  void startup(const Eigen::MatrixXd& X, int maxN) {
    linkFunction->startup(getLP(X), maxN);
    logPost = 0.;
  }
  double considerCandidate(Eigen::VectorXd& c) {
    double lp = getLP(c)(0);
    linkFunction->setCandidate(lp, &c);
    transformedLP = (*linkFunction)(lp, true);
    return transformedLP;
  }
  void acceptCandidate(bool success) {
    linkFunction->acceptCandidate(success);
    if (success) logPost += transformedLP; else
      logPost += log1p(-exp(transformedLP));
  }
  // Returns the log-posterior contribution
  double wrapup() {
    parameters = linkFunction->wrapup(prior);
    return logPost + prior->logDensity(parameters);
  }

  // Getter
  Eigen::VectorXd getParameters() {return parameters;}

  // Constructors & destructor
  BinaryRegression() {}
  BinaryRegression(Prior* p, LinkFunction* lf, Eigen::VectorXd init) :
    prior(p), linkFunction(lf), parameters(init) {}
  virtual ~BinaryRegression() {
    delete prior;
    delete linkFunction;
  };
};

class LinearRegression : public BinaryRegression {
  Eigen::VectorXd getLP(const Eigen::MatrixXd& covariates) {
    return covariates * parameters;
  }
public:
  LinearRegression(Prior* p, LinkFunction* lf, Eigen::VectorXd init) :
  BinaryRegression(p, lf, init) {}
};

#endif
