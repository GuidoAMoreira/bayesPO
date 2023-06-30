#include <RcppEigen.h>
#include "full_conditionals.h"
#include "PolyaGamma.h"
#include "utils.h"

// Constructors
gamma_prior::gamma_prior(Rcpp::List pars)
{
  shape = pars["a"];
  rate = pars["b"];
}

logit_normal::logit_normal(Rcpp::List pars, Eigen::MatrixXd x) :
  BetaDelta(x)
{
  mu = pars["mean"];
  Sigma = pars["covariance"]; // Assumed precision matrix
  Bb = Sigma * mu;
  s = mu.size();
}

BetaDelta::BetaDelta(Eigen::MatrixXd x) :
  xPositions(x) {}

// updaters
double gamma_prior::update(long ns, double area)
{
  double out, a = shape + ns, b = rate + area;
  l = R::rgamma(a, 1 / b);

  // log posterior contribution
  out = -l * b + std::log(l) * (a - 1);

  return out;
}

// Algorithm described in Polson et al. (2012)
double logit_normal::update(const Eigen::MatrixXd& onesCov,
                            const Eigen::MatrixXd& zerosCov)
{
  unsigned long i, n1 = onesCov.rows(), n0 = zerosCov.rows();
  PolyaGamma pg(1);
  Eigen::MatrixXd V = Sigma, x1 = Eigen::MatrixXd(n1, s),
    x0 = Eigen::MatrixXd(n0, s);
  Eigen::VectorXd med = Bb, xb1(n1), xb0(n0), emm;
  x1.leftCols(1) = Eigen::MatrixXd::Constant(n1,1,1);
  x1.rightCols(s - 1) = onesCov;
  x0.leftCols(1) = Eigen::MatrixXd::Constant(n0,1,1);
  x0.rightCols(s - 1) = zerosCov;
  xb1 = x1 * effects;
  xb0 = x0 * effects;

  // Calculating X' Omega X + B and X' kappa + B b
  for (i = 0; i < n1; i++) // From the data matrix X
  {
    V += pg.draw_like_devroye(xb1(i)) * x1.row(i).transpose() * x1.row(i);
    med += x1.row(i) * 0.5;
  }
  for (i = 0; i < n0; i++) // From the data matrix X
  {
    V += pg.draw_like_devroye(xb0(i)) * x0.row(i).transpose() * x0.row(i);
    med -= x0.row(i) * 0.5;
  }
  V = V.inverse();
  effects = sampleNormal(V) + V * med;

  // log posterior contribution
  emm = effects - mu;
  double out = -0.5 * (emm.transpose() * Sigma * emm)(0);
  out += link(x1, false).sum();
  out += link(x0, true).sum();
  return out;
}

// Other methods
Eigen::MatrixXd logit_normal::link(const Eigen::MatrixXd& covs, bool complementary)
{
  Eigen::MatrixXd oneCov(covs.rows(), s);


  if (covs.cols() < s)
  {
    oneCov.leftCols(1) = Eigen::MatrixXd::Constant(covs.rows(), 1, 1);
    oneCov.rightCols(covs.cols()) = covs;
  }
  else
    oneCov = covs;

  if (complementary)
    return -(oneCov * effects).array().exp().log1p();
  else
    return -(-(oneCov * effects)).array().exp().log1p();
}
