#include "markovchain.h"
extern "C" {
#include "safeR.h"
}
#include "samplers.h"
#include <omp.h>
#include "covariates.h"

// Update Markov Chain
void MarkovChain::update()
{
  applyTransitionKernel();
  iteration++;
}

void MarkovChain::update(int times)
  {for (int i=0; i < times; i++) {update();}}

void MarkovChain::applyTransitionKernel()
{
  GetRNGstate();
  logPosterior = fullConditionals() + updateLambdaStar();
  PutRNGstate();
}

double MarkovChain::updateLambdaStar()
{
  double a = lambdaA + U.size() + Xprime.size() + X.size(), b = lambdaB + area;
  lambdaStar = R::rgamma(a, 1 / b);

  // log posterior contribution
  return -lambdaStar * b + std::log(lambdaStar) * (a - 1);
}

// Latent processes full conditional
double MarkovChain::fullConditionals()
{
  // Selecting candidate cells and associating them with processes
  int nTotal = R::rpois(area * lambdaStar), iXp = 0, iU = 0, candidate;
  if (!nTotal) {
    U.resize(0);
    Xprime.resize(0);
    return 0.;
  }

  U.resize(nTotal);
  Xprime.resize(nTotal);
  Eigen::VectorXd candidateInt, candidateObs;
  double q, unif;

  intensityRegression->startup(zX, nTotal);
  observabilityRegression->startup(wX, nTotal);

#pragma omp parallel for private(candidate, candidateInt, q, unif)
  for (int i = 0; i < nTotal; i++) {
    candidate = background->pickRandomPoint();
    candidateInt = background->retrieveInt(candidate);
    unif = std::log(runif(0., 1.));
    q = intensityRegression->considerCandidate(candidateInt);
    if (unif > q) { // Accept U
#pragma omp critical
      U[iU++] = candidate;
      intensityRegression->acceptCandidate(false);
    } else {
      candidateObs = background->retrieveObs(candidate); // Stored because pointer is passed inside
      if (unif > q + observabilityRegression->considerCandidate(candidateObs)) { // Accept Xprime
#pragma omp critical
        Xprime[iXp++] = candidate;
        intensityRegression->acceptCandidate(true);
        observabilityRegression->acceptCandidate(false);
        background->addAcceptedXprime(candidate); // Add 1 to the respective cell counter. This is for the heat map.
      }
    }
  }

  U.resize(iU);
  Xprime.resize(iXp);

  return - lgamma(iXp + 1) - lgamma(iU + 1) +
    intensityRegression->wrapup() + observabilityRegression->wrapup();
}
