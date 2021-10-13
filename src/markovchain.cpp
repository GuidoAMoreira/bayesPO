#include "markovchain.h"

// Constructor
mcStep::mcStep(Eigen::VectorXd b, Eigen::VectorXd d, double l, retrievCovs *bb,
               double a, std::vector<long> x, Eigen::MatrixXd zx,
               Eigen::MatrixXd wx) :
  area(a), X(x), zX(zx), wX(wx),iteration(1), background(bb),
  ibeta(b), idelta(d), ilambdaStar(l) {}


// Update Markov Chain
void mcStep::update()
{
  applyTransitionKernel();
  iteration++;
}

void mcStep::update(long times)
{for (long i=0; i < times; i++) {update();}}

void mcStep::applyTransitionKernel()
{
  GetRNGstate();
  logPosterior = FullConditional_processes();
  logPosterior += lambda->update(X.size() + Xprime.size() + U.size(), area);
  logPosterior += beta->update(zXXp, zU);
  logPosterior += delta->update(wX, wXp);
  PutRNGstate();
}

// Latent processes full conditional
double mcStep::FullConditional_processes()
{
  // Selecting candidate cells and associating them with processes
  long nTotal = R::rpois(area * lambda->l), i, candidate, iXp = 0, iU = 0;
  std::vector<long> temp(nTotal);
  double unif, q, p;
  Eigen::VectorXd tempInt(beta->s - 1), tempObs(delta->s - 1);

  for (i = 0; i < nTotal; i++){
    candidate = background->pickRandomPoint();
    // If the uniform is larger than q, then the point is U.
    // Else, if it is larger than pq, the point is X'.
    unif = std::log(R::runif(0, 1));
    tempInt = background->retrieveInt(candidate);
    q = beta->link(tempInt.transpose(), false)(0);
    if (unif > q){ // Accept U
      temp[iU++] = candidate;
    } else { // Reject U
      tempObs = background->retrieveObs(candidate);
      p = delta->link(tempObs.transpose(), false)(0);
      if (unif > p + q) { // Accept X'
        temp[nTotal - 1 - iXp++] = candidate;
      }
    } // If neither U nor X' is accepted, point is discarded
  }

  U = std::vector<long>(&temp[0], &temp[iU]);
  Xprime = std::vector<long>(&temp[nTotal - iXp], &temp[nTotal]);
  zXXp.resize(X.size() + iXp, beta->s - 1);
  zXXp.topRows(X.size()) = zX;
  if (iXp){
    background->putInt(zXXp, temp, nTotal - iXp, nTotal - 1, X.size());
    wXp.resize(iXp, delta->s - 1);
    background->putObs(wXp, temp, nTotal - iXp, nTotal - 1, 0);
  }
  else
    wXp = Eigen::MatrixXd(0, delta->s - 1);
  if (iU){
    zU.resize(iU, beta->s - 1);
    background->putInt(zU, temp, 0, iU - 1, 0);
  }
  else
    zU = Eigen::MatrixXd(0, beta->s - 1);

  return - lgamma(iXp + 1) - lgamma(iU + 1);
}
