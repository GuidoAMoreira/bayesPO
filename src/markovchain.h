#ifndef __MARKOVCHAIN_BAYESPO_H__
#define __MARKOVCHAIN_BAYESPO_H__

#include <RcppEigen.h>
#include "covariates.h"
#include "binary_regression.h"

class MarkovChain {
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
  const double lambdaA, lambdaB;
  double lambdaStar;
  BinaryRegression* intensityRegression;
  BinaryRegression* observabilityRegression;

  double updateLambdaStar();

  void applyTransitionKernel();
  virtual double fullConditionals();
public:
  // Constructor & destructor
  MarkovChain(double a, std::vector<int> x, Eigen::MatrixXd& zx,
              Eigen::MatrixXd& wx, retrievCovs* bg, double la, double lb,
              double lambda, BinaryRegression* Int, BinaryRegression* Obs) :
  area(a), X(x), zX(zx), wX(wx), iteration(0), background(bg), lambdaA(la),
  lambdaB(lb), lambdaStar(lambda), intensityRegression(Int),
  observabilityRegression(Obs) {}
  ~MarkovChain() {
    delete background;
    delete intensityRegression;
    delete observabilityRegression;
  }

  // Methods
  void update();
  void update(int times);

  // Getters
  Eigen::VectorXd getBeta() {return intensityRegression->getParameters();}
  Eigen::VectorXd getDelta() {return observabilityRegression->getParameters();}
  double getLambda() {return lambdaStar;}
  double getLogPosterior() {return logPosterior;}
  int getUsize() {return U.size();}
  int getXprimeSize() {return Xprime.size();}
  Eigen::VectorXd getHeatMap() {return background->getUnobservedCounts();}
};


#endif
