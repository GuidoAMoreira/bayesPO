#include "markovchain.h"
#include "safeR.h"

// Constructor
mcStep::mcStep(Eigen::VectorXd b, Eigen::VectorXd d, double l, retrievCovs *bb,
               double a, std::vector<int> x, Eigen::MatrixXd zx,
               Eigen::MatrixXd wx) :
  area(a), X(x), zX(zx), wX(wx), zXXp(zx), iteration(1), background(bb),
  ibeta(b), idelta(d), ilambdaStar(l) {}


// Update Markov Chain
void mcStep::update()
{
  applyTransitionKernel();
  iteration++;
}

void mcStep::update(int times)
{for (int i=0; i < times; i++) {update();}}

void mcStep::applyTransitionKernel()
{
  GetRNGstate();
  logPosterior = FullConditional_processes();
  logPosterior += lambda->update(X.size() + Xprime.size() + U.size(), area);
  logPosterior += beta->update(zXXp, zU);
  logPosterior += delta->update(wX, wXp);
  PutRNGstate();
}

// // Latent processes full conditional
// double mcStep::FullConditional_processes()
// {
//   // Selecting candidate cells and associating them with processes
//   int nTotal = R::rpois(area * lambda->l), i, candidate, iXp = 0, iU = 0;
//   std::vector<int> temp(nTotal);
//   double unif, q, p;
//   Eigen::VectorXd tempInt(beta->s - 1), tempObs(delta->s - 1);
//
//   for (i = 0; i < nTotal; i++){
//     candidate = background->pickRandomPoint();
//     // If the uniform is larger than q, then the point is U.
//     // Else, if it is larger than pq, the point is X'.
//     unif = std::log(R::runif(0, 1));
//     tempInt = background->retrieveInt(candidate);
//     q = beta->link(tempInt.transpose(), false)(0);
//     if (unif > q){ // Accept U
//       temp[iU++] = candidate;
//     } else { // Reject U
//       tempObs = background->retrieveObs(candidate);
//       p = delta->link(tempObs.transpose(), false)(0);
//       if (unif > p + q) { // Accept X'
//         temp[nTotal - 1 - iXp++] = candidate;
//         background->addAcceptedXprime(candidate);
//       }
//     } // If neither U nor X' is accepted, point is discarded
//   }
//
//   U = std::vector<int>(&temp[0], &temp[iU]);
//   Xprime = std::vector<int>(&temp[nTotal - iXp], &temp[nTotal]);
//   zXXp.resize(X.size() + iXp, beta->s - 1);
//   zXXp.topRows(X.size()) = zX;
//   if (iXp){
//     background->putInt(zXXp, temp, nTotal - iXp, nTotal - 1, X.size());
//     wXp.resize(iXp, delta->s - 1);
//     background->putObs(wXp, temp, nTotal - iXp, nTotal - 1, 0);
//   }
//   else
//     wXp = Eigen::MatrixXd(0, delta->s - 1);
//   if (iU){
//     zU.resize(iU, beta->s - 1);
//     background->putInt(zU, temp, 0, iU - 1, 0);
//   }
//   else
//     zU = Eigen::MatrixXd(0, beta->s - 1);
//
//   return - lgamma(iXp + 1) - lgamma(iU + 1);
// }

// Latent processes full conditional
double mcStep::FullConditional_processes()
{
  // Selecting candidate cells and associating them with processes
  int nTotal = R::rpois(area * lambda->l), i, notU = 0, notXp = 0, iXp = 0, iU = 0;//, candidate;
  //std::vector<int> temp(nTotal);
  //double unif, q, p;
  //Eigen::VectorXd tempInt(beta->s - 1), tempObs(delta->s - 1);
  if (!nTotal) {
    zU = Eigen::MatrixXd(0, beta->s - 1);
    zXXp.conservativeResize(zX.rows() + iXp, beta->s - 1);
    wXp = Eigen::MatrixXd(0, delta->s - 1);
    return 0;
  }
  Eigen::VectorXi candidates = background->pickRandomPoint(nTotal);
  Eigen::MatrixXd candidateInt = background->retrieveInt(candidates);

  // U process
  Eigen::VectorXd selectors = runif(nTotal).array().log().matrix() -
    beta->link(candidateInt, false);
  Eigen::VectorXi permIdx(nTotal);
  for (i = 0; i < nTotal; i++)
    if (selectors(i) > 0) permIdx(iU++) = i; else permIdx(nTotal - ++notU) = i;

  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
  perm.indices() = permIdx;
  perm = perm.transpose().eval();
  candidateInt = perm * candidateInt;
  candidates = perm * candidates;
  selectors = (perm * selectors).tail(notU);
  if (iU) zU = candidateInt.topRows(iU); else
    zU = Eigen::MatrixXd(0, beta->s - 1);
  U = std::vector<int>(candidates.data(), candidates.data() + iU);

  // Xp process
  nTotal -= iU;
  Eigen::VectorXi candidates2 = candidates.tail(notU);
  permIdx.resize(nTotal);
  Eigen::MatrixXd candidateObs = background->retrieveObs(candidates2);
  selectors -= delta->link(candidateObs, false);
  for (i = 0; i < nTotal; i++)
    if (selectors(i) > 0) permIdx(iXp++) = i; else
      permIdx(nTotal - ++notXp) = i;

  zXXp.conservativeResize(zX.rows() + iXp, beta->s - 1);
  if (iXp) {
    perm.indices() = permIdx;
    perm = perm.transpose().eval();
    Eigen::MatrixXd temp = (perm * candidateInt.bottomRows(notU)).topRows(iXp);
    zXXp.bottomRows(iXp) = temp;
    wXp = (perm * candidateObs).topRows(iXp);
    candidates = perm * candidates2;
    background->addAcceptedXprime(candidates.head(iXp));
  } else wXp = Eigen::MatrixXd(0, delta->s - 1);
  Xprime = std::vector<int>(candidates.data(), candidates.data() + iXp);

  /*
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
        background->addAcceptedXprime(candidate);
      }
    } // If neither U nor X' is accepted, point is discarded
  }

  U = std::vector<int>(&temp[0], &temp[iU]);
  Xprime = std::vector<int>(&temp[nTotal - iXp], &temp[nTotal]);
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
   */

  return - lgamma(iXp + 1) - lgamma(iU + 1);
}
