#include "markovchain.h"

// Constructor
mcStep::mcStep(NumericVector b, NumericVector d, double l, retrievCovs *bb,
               double a, std::vector<R_xlen_t> x, mat zx, mat wx) :
  nb(b.length()),nd(d.length()),area(a),X(x),zX(zx),beta(vec(b)),delta(vec(d)),
  lambdaStar(l),wX(wx),iteration(1),background(bb),ibeta(vec(b)),idelta(vec(d)),
  ilambdaStar(l) {}

// Update Markov Chain
void mcStep::update()
{
  applyTransitionKernel();
  iteration++;
}

// Attaching function depending on model
void mcStep::attachLambdaFC(double FCl(double *l, R_xlen_t x, R_xlen_t xp,
                                       R_xlen_t u, List p, double a))
{
  FullConditional_lambda = FCl;
}

//void mcStep::attachBetaFC(double FCb(vec *b, mat *oc, mat *zc, List p, double link(arma::vec b, arma::vec cov, bool c)))
void mcStep::attachBetaFC(double FCb(vec *b, mat *oc, mat *zc, List p,
                                     double link(vec b, vec cov, bool c)))
{
  FullConditional_beta = FCb;
}

//void mcStep::attachDeltaFC(double FCd(vec *d, mat *oc, mat *zc, List p, double link(arma::vec b, arma::vec cov, bool c)))
void mcStep::attachDeltaFC(double FCd(vec *d, mat *oc, mat *zc, List p,
                                      double link(vec b, vec cov, bool c)))
{
  FullConditional_delta = FCd;
}

void mcStep::attachCalcq(double cq(vec b, vec co, bool c))
{
  calc_q = cq;
}

void mcStep::attachCalcp(double cp(vec d, vec co, bool c))
{
  calc_p = cp;
}

void mcStep::applyTransitionKernel()
{
  GetRNGstate();
  logPosterior = FullConditional_processes();
  logPosterior += (*FullConditional_lambda)(&lambdaStar, X.size(), Xprime.size(),
                   U.size(), lambdaPars, area);
  logPosterior += (*FullConditional_beta)(&beta, &zXXp, &zU, betaPars, calc_q);
  logPosterior += (*FullConditional_delta)(&delta, &wX, &wXp, deltaPars, calc_p);
  PutRNGstate();
}

// Latent processes full conditional
double mcStep::FullConditional_processes()
{
  // Selecting candidate cells and associating them with processes
  R_xlen_t nTotal = R::rpois(area*lambdaStar), i, candidate, iXp = 0, iU = 0;
  Xprime = std::vector<R_xlen_t>();
  U = std::vector<R_xlen_t>();
  double unif, q, p;
  vec tempInt(nb-1), tempObs(nd-1);
  mat temp_wXp(nTotal,nd-1), temp_zXp(nTotal,nb-1), temp_zU(nTotal,nb-1);

  for (i=0;i<nTotal;i++){
    candidate = R_xlen_t(R::runif(0,1) * (*background).ncell);
    // If the uniform is larger than q, then the point is U.
    // Else, if it is larger than pq, the point is X'.
    unif = std::log(R::runif(0,1));
    tempInt = (*background).retrieveInt(candidate);
    q = calc_q(beta,tempInt,false);
    if (unif > q){ // Accept U
      U.push_back(candidate);
      temp_zU.row(iU++) = tempInt.t();
    } else { // Reject U
      tempObs = (*background).retrieveObs(candidate);
      p = calc_p(delta,tempObs,false);
      if (unif > p + q) { // Accept X'
        Xprime.push_back(candidate);
        temp_zXp.row(iXp) = tempInt.t();
        temp_wXp.row(iXp++) = tempObs.t();
      }
    } // If neither U nor X' is accepted, point is discarded
  }

  zXXp = join_vert(zX,temp_zXp.head_rows(Xprime.size()));
  wXp = temp_wXp.head_rows(Xprime.size());
  zU = temp_zU.head_rows(U.size());

  return  - lgamma(Xprime.size()) - lgamma(U.size());
}

