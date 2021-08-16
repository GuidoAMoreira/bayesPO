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

logit_normal::logit_normal(Rcpp::List pars, Eigen::MatrixXd x, Covariance *c) :
  BetaDelta(x, c)
{
  mu = pars["mean"];
  Sigma = pars["covariance"]; // Assumed precision matrix
  Bb = Sigma * mu;
  s = mu.size();
  xMarks = Eigen::ArrayXd::Zero(xPositions.rows()).matrix();
}

BetaDelta::BetaDelta(Eigen::MatrixXd x, Covariance *c) :
  xPositions(x), covariance(c) {
  long n = xPositions.rows(), i, j;
  xDists = Eigen::MatrixXd::Identity(n, n);
  double dist;

  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
    {
      dist = sqrt((xPositions.row(i) * xPositions.row(i).transpose())(0));
      xDists(i, j) = dist;
      xDists(j, i) = dist;
    }
}

Covariance::Covariance(bool use) : used(use) {}

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
  Eigen::VectorXd med = mu, xb1(n1), xb0(n0), emm;
  x1.leftCols(1) = Eigen::MatrixXd::Constant(n1,1,1);
  x1.rightCols(s - 1) = onesCov;
  x0.leftCols(1) = Eigen::MatrixXd::Constant(n0,1,1);
  x0.rightCols(s - 1) = zerosCov;
  xb1 = x1 * effects;
  xb0 = x0 * effects;

  // Calculating X' Omega X + B and X' kappa + B b
  // offset represents a Gaussian Process if it is in the model
  for (i = 0; i < n1; i++) // From the data matrix X
  {
    //polyagamma[i] = pg.draw_like_devroye(xb1(i) + offset(i));
    V += pg.draw_like_devroye(xb1(i)) * x1.row(i).transpose() * x1.row(i);
    //med += x1.row(i) * (0.5 - offset(i) * polyagamma(i));
    med += x1.row(i) * 0.5;
  }
  for (i = 0; i < n0; i++) // From the data matrix X
  {
    //polyagamma[i + n1] = pg.draw_like_devroye(xb0(i) + offset(i + n1));
    V += pg.draw_like_devroye(xb0(i)) * x0.row(i).transpose() * x0.row(i);
    //med -= x0.row(i) * (0.5 + offset(i + n1) * polyagamma(i + n1));
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

double BetaDelta::updateCovariance()
{return covariance->update(this);}

// INCOMPLETE
double logit_normal::updateGPx(const Eigen::MatrixXd& xCovar)
{
  if (!covariance->used) return 0.;

  Eigen::MatrixXd covars;
  if (initializedGP)
  {
    Eigen::VectorXd muX = polyagamma.head(xPositions.rows()).array() *
      (xCovar * effects).array() - 0.5;

    return -0.5 * xMarks.transpose() * covars.inverse() * xMarks;
  }

  // First iteration: Just sample from the prior
  covars = covariance->applyCovariance(xDists);
  xMarks = sampleNormal(covars);
  initializedGP = true;
  return -0.5 * xMarks.transpose() * covars.inverse() * xMarks;
}

// INCOMPLETE
double logit_normal::updateGP(const Eigen::MatrixXd &pos,
                              const Eigen::MatrixXd &onesC,
                              const Eigen::MatrixXd &zerosC)
{
  if (!covariance->used) // Skip if Gaussian Process is not used in the model
  {
    offset = Eigen::ArrayXd::Zero(onesC.rows() + zerosC.rows()).matrix();
    return 0.;
  }

  Positions = pos;
  if (initializedGP)
  {

  }

  // First iteration

  return 0.;
}

// INCOMPLETE
double Exponential_Covariance::update(BetaDelta* points)
{
  double proposal = R::rnorm(phi, 0.3), temp = phi, propDens, pastDens, alpha;
  pastDens = - rate * phi + (shape - 1) * log(phi);
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
    //return (oneCov * effects + offset.bottomRows(covs.rows())).array().exp().log1p();
    return -(oneCov * effects).array().exp().log1p();
  else
    //return (-(oneCov * effects + offset.topRows(covs.rows()))).array().exp().log1p();
    return -(-(oneCov * effects)).array().exp().log1p();
}

const Eigen::MatrixXd& Covariance::calcCovMatrix(const Eigen::MatrixXd& x,
                                                 const Eigen::MatrixXd& y)
{
  long nx = x.rows(), ny = y.rows(), i, j;
  Eigen::MatrixXd points(nx * ny, 2);

  for (i = 0; i < nx; i++)
    for (j = 0; i < ny; j++)
      points.row(i * ny + j) = x.row(i) - y.row(j);

  return calcCovMatrixPoints(points);
}

const Eigen::MatrixXd& Covariance::calcCovMatrixPoints(const Eigen::MatrixXd& points)
{
  R_xlen_t i, j;
  long n = points.rows();
  Eigen::MatrixXd output = Eigen::MatrixXd::Identity(n, n);
  double cov;

  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
    {
      cov = calcCov(sqrt((points.row(i) * points.row(i).transpose())(0)));
      output(i, j) = cov;
      output(j, i) = cov;
    }

  return output;
}

const Eigen::MatrixXd& Covariance::applyCovariance(const Eigen::MatrixXd& mat)
{
  long n = mat.rows();
  Eigen::MatrixXd output = Eigen::MatrixXd::Identity(n, n);
  double cov;

  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
    {
      cov = calcCov(mat(i, j));
      output(i, j) = cov;
      output(j, i) = cov;
    }

  return output;
}

double Exponential_Covariance::calcCov(double dist)
{return exp(- dist / phi) * sigma2;}

double Gaussian_Covariance::calcCov(double dist)
{return exp(- dist * dist / phi) * sigma2;}

double Matern_Covariance::calcCov(double dist)
{
  double k = LOG_SQRT_2 + log(sqrt(nu)) + log(dist) - log(rho);
  return exp(log(sigma2) + (1 - nu) * M_LN2 - lgamma(nu) +
             nu * k + log(R::bessel_k(exp(k), nu, 1)));
}
