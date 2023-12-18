#ifndef __MARKOVCHAIN_BAYESPO_H__
#define __MARKOVCHAIN_BAYESPO_H__

#include <RcppEigen.h>
#include "covariates.h"
#include "full_conditionals.h"

// Generic Markov Chain step class
// Parametric form for q(.) and p(.) are attached before MC procedure begins
class mcStep
{
public:
  //// Unchanging attributes
  const double area;
  const std::vector<int> X;
  const Eigen::MatrixXd zX;

  //// Attributes
  std::vector<int> U, Xprime; // vectors of indexes
  double logPosterior;
  //mat zX, zU, wX, wXp, zXXp;
  Eigen::MatrixXd zXXp, wXp, zU, wX;
  unsigned int iteration;
  retrievCovs *background;

  // Initialization
  Eigen::VectorXd ibeta, idelta;
  double ilambdaStar;

  // Parameters and their updaters
  LambdaStar *lambda;
  BetaDelta *beta, *delta;

  // Constructor
  mcStep(Eigen::VectorXd b, Eigen::VectorXd d, double l, retrievCovs *bb,
         double a, std::vector<int> x, Eigen::MatrixXd zx,
         Eigen::MatrixXd wx);

  // Methods
  void update();
  void update(int times);

protected:
  void applyTransitionKernel(); // Markov chain transition kernel

  // Full conditionals
  double FullConditional_processes();
};

/*
// Class to use a Gaussian process in the intensity set
class mcStep_GP_int : public mcStep
{
public:
  Eigen::VectorXd Zi;

  mcStep_GP_int(Eigen::VectorXd b, Eigen::VectorXd d, double l, retrievCovs *bb,
                double a, std::vector<int> x, Eigen::MatrixXd zx,
                Eigen::MatrixXd wx);

private:
  void applyTransitionKernel(); // Markov chain transition kernel

  // Full conditionals
  double FullConditional_processes();
};
*/
#endif
