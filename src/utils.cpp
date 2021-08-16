#include "utils.h"
#include "full_conditionals.h"

mcStep startup_mcmc(Eigen::VectorXd beta, Eigen::VectorXd delta, double lambda,
                    std::string b_updater, std::string d_updater,
                    std::string l_updater, Rcpp::List parB, Rcpp::List parD,
                    Rcpp::List parL, std::string xClass,
                    Eigen::MatrixXd xValues, Eigen::VectorXi xI,
                    Eigen::VectorXi xO, retrievCovs *covars, double aD)
{

  // Data structures
  std::vector<long> x;
  Eigen::MatrixXd zx, wx;

  // Data matrix
  if (xClass == "num_mat") {Eigen::MatrixXd X_covs = xValues;
    importX_double(X_covs, beta.size(), delta.size(), xI, xO, x, zx, wx);}
  else if (xClass == "reference") {Eigen::VectorXd X_covs = xValues;
    determineX(X_covs, beta.size(), delta.size(), covars, x, zx, wx);}

  // Creating and initializing the Markov Chain
  mcStep mc(beta, delta, lambda, covars, aD, x, zx, wx);

  // Attaching Lambda updater
  if (l_updater == "gamma")
      mc.lambda = new gamma_prior(parL);
  mc.lambda->l = lambda;

  // Attaching Beta updater
  if (b_updater == "logit_normal")
    mc.beta = new logit_normal(parB, Eigen::MatrixXd(), new Covariance(false));
  mc.beta->effects = beta;

  // Attaching Delta updater
  if (d_updater == "logit_normal")
    mc.delta = new logit_normal(parD, Eigen::MatrixXd(), new Covariance(false));
  mc.delta->effects = delta;

  return mc;
}

// Functions to determine data X when it comes with a matrix with the same number of columns as
// the background
void importX_double(Eigen::MatrixXd x, long nb, long nd,
                    Eigen::VectorXi xI, Eigen::VectorXi xO,
                    std::vector<long> &x_data, Eigen::MatrixXd &zx_data,
                    Eigen::MatrixXd &wx_data)
{
  x_data = std::vector<long>(x.rows());
  Eigen::MatrixXd zX(x.rows(), nb - 1), wX(x.rows(), nd - 1);
  long i, j;
  for (i = 0;i<x.rows();i++)
  {
    for (j = 0; j < (nb - 1); j++)
    {
      zX(i, j) = x(i, xI[j]);
    }
    for (j = 0; j < (nd - 1); j++)
    {
      wX(i, j) = x(i, xO[j]);
    }
  }

  zx_data = zX;
  wx_data = wX;
}

void importX_int(Eigen::VectorXi x, long nb, long nd,
                 Eigen::VectorXi xI, Eigen::VectorXi xO,
                 std::vector<long> &x_data, Eigen::MatrixXd &zx_data,
                 Eigen::MatrixXd &wx_data)
{
  x_data = std::vector<long>(x.rows());
  Eigen::MatrixXd zX(x.rows(),nb-1), wX(x.rows(),nd-1);
  long i, j;
  for (i = 0; i < x.rows(); i++)
  {
    for (j = 0; j < (nb - 1); j++)
    {
      zX(i, j) = double(x(i, xI[j]));
    }
    for (j = 0; j < (nd - 1); j++)
    {
      wX(i, j) = double(x(i, xO[j]));
    }
  }

  zx_data = zX;
  wx_data = wX;
}

// Function to determine data X when it comes with a vector of indexes relative
// to the background cells
void determineX(Eigen::VectorXd x, long nb,
                      long nd, retrievCovs *b,
                      std::vector<long> &x_data, Eigen::MatrixXd &zx_data,
                      Eigen::MatrixXd &wx_data)
{
  std::vector<long> X(x.size());
  Eigen::MatrixXd zX(x.size(),nb - 1), wX(x.size(), nd - 1);
  long i;
  for (i=0; i < x.size() ; i++) {
    X.push_back(x[i]);
    zX.row(i) = (*b).retrieveInt(long(x[i]));
    wX.row(i) = (*b).retrieveObs(long(x[i]));
  }

  x_data = X;
  zx_data = zX;
  wx_data = wX;
}
