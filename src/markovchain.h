#ifndef __MARKOVCHAIN_BAYESPO_H__
#define __MARKOVCHAIN_BAYESPO_H__

#include <RcppEigen.h>
#include "covariates.h"
#include "binary_regression.h"
#include "full_conditionals.h"

// Generic Markov Chain step class
// Parametric form for q(.) and p(.) are attached before MC procedure begins
class MarkovChain
{
  //// Unchanging attributes
  const double area;
  const std::vector<int> X;
  const Eigen::MatrixXd zX;
  const Eigen::MatrixXd wX;

  //// Process variables
  std::vector<int> U, Xprime; // vectors of indexes
  double logPosterior;
  unsigned int iteration;
  retrievCovs *background;

  // Parameters
  Eigen::VectorXd beta, delta;
  double lambda;

  // Prior
  RegressionPrior* betaPrior, deltaPrior;

  void applyTransitionKernel(); // Markov chain transition kernel
public:
  // Constructor
  MarkovChain(Eigen::VectorXd b, Eigen::VectorXd d, double l, retrievCovs *bb,
              double a, std::vector<int>& x, Eigen::MatrixXd& zx,
              Eigen::MatrixXd& wx) : area(a), X(x), zX(zx), wX(wx),
              iteration(0), background(bb), beta(b), delta(d), lambda(l) {}

  // Methods
  void update();
  void update(int times);

protected:
  virtual double FullConditionals();
};


class LogisticLogistic : public MarkovChain {
  double FullConditionals();
};

#endif
