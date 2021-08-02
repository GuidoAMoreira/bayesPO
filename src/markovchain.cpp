#include "markovchain.h"

// Constructor
mcStep::mcStep(Eigen::VectorXd b, Eigen::VectorXd d, double l, retrievCovs *bb,
               double a, std::vector<long> x, Eigen::MatrixXd zx,
               Eigen::MatrixXd wx) :
  area(a), X(x), zX(zx), wX(wx),iteration(1), background(bb),
  ibeta(b), idelta(d), ilambdaStar(l) {}

/*mcStep_GP_int::mcStep_GP_int(Eigen::VectorXd b, Eigen::VectorXd d, double l,
                             retrievCovs *bb, double a,
                             std::vector<long> x, Eigen::MatrixXd zx,
                             Eigen::MatrixXd wx) :
  mcStep(b, d, l, bb, a, x, zx, wx) {}*/

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
  logPosterior += delta->update(wXp, wXp);
  PutRNGstate();
}

// Latent processes full conditional
double mcStep::FullConditional_processes()
{
  // Selecting candidate cells and associating them with processes
  long nTotal = R::rpois(area * lambda->l), i, candidate, iXp = 0, iU = 0;
  Xprime = std::vector<long>();
  U = std::vector<long>();
  std::vector<long> Utemp(nTotal), Xptemp(nTotal);
  double unif, q, p;
  Eigen::VectorXd tempInt(beta->s - 1), tempObs(delta->s - 1);
  Eigen::MatrixXd temp_wXp(nTotal, delta->s - 1), temp_zXp(nTotal, beta->s - 1),
      temp_zU(nTotal, beta->s - 1);

  for (i = 0; i < nTotal; i++){
    candidate = background->pickRandomPoint();
    // If the uniform is larger than q, then the point is U.
    // Else, if it is larger than pq, the point is X'.
    unif = std::log(R::runif(0, 1));
    tempInt = background->retrieveInt(candidate);
    q = beta->link(tempInt.transpose(), false)(0);
    if (unif > q){ // Accept U
      Utemp[iU] = candidate;
      temp_zU.row(iU++) = tempInt;
    } else { // Reject U
      tempObs = background->retrieveObs(candidate);
      p = delta->link(tempObs.transpose(), false)(0);
      if (unif > p + q) { // Accept X'
        Xptemp[iXp] = candidate;
        temp_zXp.row(iXp) = tempInt;
        temp_wXp.row(iXp++) = tempObs;
      }
    } // If neither U nor X' is accepted, point is discarded
  }

  U = std::vector<long>(&Utemp[0], &Utemp[iU]);
  Xprime = std::vector<long>(&Xptemp[0], &Xptemp[iXp]);
  zXXp.resize(X.size() + iXp, beta->s - 1);
  zXXp.topRows(X.size()) = zX;
  zXXp.bottomRows(iXp) = temp_zXp.topRows(iXp);
  wXp.resize(iXp, delta->s - 1);
  wXp = temp_wXp.topRows(iXp);
  zU.resize(iU, beta->s - 1);
  zU = temp_zU.topRows(iU);

  return - lgamma(iXp + 1) - lgamma(iU + 1);
}
/*
// Gaussian Process in the intensity predictor
void mcStep_GP_int::applyTransitionKernel()
{
  GetRNGstate();
  logPosterior = FullConditional_processes();
  logPosterior += (*FullConditional_lambda)(lambdaStar, X.size(), Xprime.size(),
                   U.size(), lsPrior, area);
  logPosterior += (*FullConditional_beta)(beta, zXXp, zU, bPrior, Zi, calc_q);
  logPosterior += (*FullConditional_delta)(delta, wX, wXp, dPrior,
                   Eigen::ArrayXd::Zero(X.size() + Xprime.size()).matrix(),
                   calc_p);
  PutRNGstate();
}

// Latent processes full conditional
double mcStep_GP_int::FullConditional_processes()
{
  // Selecting candidate cells and associating them with processes
  long nTotal = R::rpois(area * lambdaStar), i, iXp = 0, iU = 0;
  Xprime = std::vector<long>();
  U = std::vector<long>();
  std::vector<long> Utemp(nTotal), Xptemp(nTotal);
  double unif, q, p;
  Eigen::VectorXd tempInt(nb - 1), tempObs(nd - 1);
  Eigen::VectorXi candidates = background->chooseCell(nTotal);
  Eigen::MatrixXd temp_wXp(nTotal, nd - 1), temp_zXp(nTotal, nb - 1),
  temp_zU(nTotal, nb - 1);

  for (i = 0; i < nTotal; i++){
    candidate = background->chooseCell();
    // If the uniform is larger than q, then the point is U.
    // Else, if it is larger than pq, the point is X'.
    unif = std::log(R::runif(0, 1));
    tempInt = background->retrieveInt(candidate);
    q = calc_q(beta, tempInt.transpose(), false)[0];
    if (unif > q){ // Accept U
      Utemp[iU] = candidate;
      temp_zU.row(iU++) = tempInt;
    } else { // Reject U
      tempObs = background->retrieveObs(candidate);
      p = calc_p(delta, tempObs.transpose(), false)[0];
      if (unif > p + q) { // Accept X'
        Xptemp[iXp] = candidate;
        temp_zXp.row(iXp) = tempInt;
        temp_wXp.row(iXp++) = tempObs;
      }
    } // If neither U nor X' is accepted, point is discarded
  }

  U = std::vector<long>(&Utemp[0], &Utemp[iU]);
  Xprime = std::vector<long>(&Xptemp[0], &Xptemp[iXp]);
  zXXp.resize(X.size() + iXp, nb - 1);
  zXXp.topRows(X.size()) = zX;
  zXXp.bottomRows(iXp) = temp_zXp.topRows(iXp);
  wXp.resize(iXp, nd - 1);
  wXp = temp_wXp.topRows(iXp);
  zU.resize(iU, nb - 1);
  zU = temp_zU.topRows(iU);

  return  - lgamma(iXp) - lgamma(iU);
}
*/
