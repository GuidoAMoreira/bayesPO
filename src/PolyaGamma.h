#ifndef __POLYAGAMMA__
#define __POLYAGAMMA__

//#include "simple_RNG_wrapper.h"
//#include "truncated_norm.h"
//#include "truncated_gamma.h"
//#include "inverse_gaussian.h"

#include <cmath>
#include <vector>
//#include <stdio.h>
#include <stdexcept>
#include <RcppArmadillo.h>
using namespace Rcpp;

using std::vector;

// The numerical accuracy of __PI will affect your distribution.
const double __PI = 3.141592653589793238462643383279502884197;
const double HALFPISQ = 0.5 * __PI * __PI;
const double FOURPISQ = 4 * __PI * __PI;
const double __TRUNC = 0.64;
const double __TRUNC_RECIP = 1.0 / __TRUNC;

class PolyaGamma
{

  // For sum of Gammas.
  int T;
  vector<double> bvec;

public:

  // Constructors.
  PolyaGamma(int trunc = 200);

  // Draw.
  // double draw(double n, double z, RNG& r);
  double draw(int n, double z);
  double draw_sum_of_gammas(double n, double z);
  double draw_like_devroye(double z);

  //void draw(MF x, double a, double z, RNG& r);
  //void draw(MF x, MF a, MF z, RNG& r);

  // Utility.
  void set_trunc(int trunc);

  // Helper.
  double a(int n, double x);
  double pigauss(double x, double Z);
  double mass_texpon(double Z, double fz);
  double rtigauss(double Z);

  static double jj_m1(double b, double z);
  static double jj_m2(double b, double z);
  static double pg_m1(double b, double z);
  static double pg_m2(double b, double z);

};

#endif
