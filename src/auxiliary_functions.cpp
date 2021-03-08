#include "auxiliary_functions.h"

mcStep startup_mcmc(NumericVector beta, NumericVector delta, double lambda,
                    std::string b_updater, std::string d_updater,
                    std::string l_updater, List parB, List parD, List parL,
                    std::string xClass, SEXP xValues, IntegerVector xI,
                    IntegerVector xO, retrievCovs *covars, double aD)
{

  List data;
  // Data matrix
  if (xClass == "num_mat") {NumericMatrix X_covs = xValues;data =
    importX_double(X_covs, beta.length(), delta.length(), xI, xO);}
  else if (xClass == "num_mat") {IntegerMatrix X_covs = xValues;data =
    importX_int(X_covs, beta.length(), delta.length(), xI, xO);}
  else if (xClass == "reference") {NumericVector X_covs = xValues;data =
    determineX(X_covs, beta.length(), delta.length(), covars);}

  // Creating and initializing the Markov Chain
  mcStep mc(beta,delta,lambda,covars,aD,data[0],data[1],data[2]);

  // Attaching Lambda updater
  if (l_updater == "gamma") {mc.attachLambdaFC(l_gamma);}

  // Attaching Beta updater
  if (b_updater == "logit_normal") {mc.attachBetaFC(bd_logit_normal);mc.attachCalcq(logistic);}

  // Attaching Delta updater
  if (d_updater == "logit_normal") {mc.attachDeltaFC(bd_logit_normal);mc.attachCalcp(logistic);}

  // Setting prior parameters
  mc.lambdaPars = parL;
  mc.betaPars = parB;
  mc.deltaPars = parD;

  return mc;
}

// Functions to determine data X when it comes with a matrix with the same number of columns as
// the background
List importX_double(NumericMatrix x, R_xlen_t nb, R_xlen_t nd, IntegerVector xI, IntegerVector xO)
{
  std::vector<R_xlen_t> X(x.nrow());
  mat zX(x.nrow(),nb-1), wX(x.nrow(),nd-1);
  R_xlen_t i,j;
  for (i = 0;i<x.nrow();i++)
  {
    for (j = 0; j<(nb-1);j++)
    {
      zX(i,j) = x(i,xI[j]);
    }
    for (j = 0; j<(nd-1);j++)
    {
      wX(i,j) = x(i,xO[j]);
    }
  }
  return List::create(X,zX,wX);
}

List importX_int(IntegerMatrix x, R_xlen_t nb, R_xlen_t nd, IntegerVector xI, IntegerVector xO)
{
  std::vector<R_xlen_t> X(x.nrow());
  mat zX(x.nrow(),nb-1), wX(x.nrow(),nd-1);
  R_xlen_t i,j;
  for (i = 0;i<x.nrow();i++)
  {
    for (j = 0; j<(nb-1);j++)
    {
      zX(i,j) = double(x(i,xI[j]));
    }
    for (j = 0; j<(nd-1);j++)
    {
      wX(i,j) = double(x(i,xO[j]));
    }
  }
  return List::create(X,zX,wX);
}

// Function to determine data X when it comes with a vector of indexes relative
// to the background cells
List determineX(NumericVector x, R_xlen_t nb, R_xlen_t nd, retrievCovs *b)
{
  std::vector<R_xlen_t> X(x.length());
  mat zX(x.length(),nb-1), wX(x.length(),nd-1);
  R_xlen_t i;
  for (i=0;i<x.length();i++) {
    X.push_back(x[i]);
    zX.row(i) = (*b).retrieveInt(R_xlen_t(x[i]));
    wX.row(i) = (*b).retrieveObs(R_xlen_t(x[i]));
  }
  return List::create(X,zX,wX);
}
