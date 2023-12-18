#ifndef __SAFER_HPP__
#define __SAFER_HPP__

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

inline double runif(double a = 0, double b = 1) {
  double output;
#pragma omp critical
  output = R::runif(a, b);
  return output;
}

inline Eigen::VectorXd runif(int n, double a = 0, double b = 1) {
  Eigen::VectorXd output(n);
#pragma omp critical
  output = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rcpp::runif(n, a, b));
  return output;
}

inline double rnorm(double m, double sd) {
  double output;
#pragma omp critical
  output = R::rnorm(m, sd);
  return output;
}

inline double rgamma(double shape, double scale) {
  double output;
#pragma omp critical
  output = R::rgamma(shape, scale);
  return output;
}

inline double rpois(double rate) {
  double output;
#pragma omp critical
  output = R::rpois(rate);
  return output;
}

inline Eigen::VectorXd rnorm(int n, double m = 0, double sd = 1) {
  Eigen::VectorXd output(n);
#pragma omp critical
  output = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rcpp::rnorm(n, m, sd));
  return output;
}

#endif
