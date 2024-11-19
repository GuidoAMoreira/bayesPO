#ifndef __UTILS_BAYESPO_H__
#define __UTILS_BAYESPO_H__

#include "covariates.h"
#include "markovchain.h"
#include <RcppEigen.h>
extern "C" {
#include "safeR.h"
}

// Formating x functions
void importX(Eigen::MatrixXd x, int nb, int nd,
                    Eigen::VectorXi xI, Eigen::VectorXi xO,
                    std::vector<int> &x_data, Eigen::MatrixXd &zx_data,
                    Eigen::MatrixXd &wx_data);

// Sampling. Simple inline function
inline Eigen::VectorXd sampleNormal(const Eigen::MatrixXd& Sigma)
{return Sigma.llt().matrixL() *
  Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rcpp::rnorm(Sigma.rows(), 0, 1));}

inline Eigen::VectorXd safe_pnorm(const Eigen::VectorXd& x,
                                  const Eigen::VectorXd& mu,
                                  const Eigen::VectorXd& sd,
                                  bool lower, bool log) {
  Eigen::VectorXd output(x.size());
  for (int i = 0; i < x.size(); i++) {
    output(i) = safe_pnorm(x(i), mu(i), sd(i), lower, log);
  }
  return output;
}

inline Eigen::VectorXd safe_pnorm(const Eigen::VectorXd& x,
                                  const Eigen::VectorXd& mu,
                                  double sd, bool lower, bool log) {
  Eigen::VectorXd output(x.size());
  for (int i = 0; i < x.size(); i++) {
    output(i) = safe_pnorm(x(i), mu(i), sd, lower, log);
  }
  return output;
}

inline Eigen::VectorXd safe_pnorm(const Eigen::VectorXd& x, double mu,
                                  double sd, bool lower, bool log) {
  Eigen::VectorXd output(x.size());
  for (int i = 0; i < x.size(); i++) {
    output(i) = safe_pnorm(x(i), mu, sd, lower, log);
  }
  return output;
}

#endif

