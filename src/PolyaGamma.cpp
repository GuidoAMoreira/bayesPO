#include "PolyaGamma.h"

////////////////////////////////////////////////////////////////////////////////
// Constructors //
////////////////////////////////////////////////////////////////////////////////

PolyaGamma::PolyaGamma(int trunc) : T(trunc), bvec(T)
{
  set_trunc(T);
} // PolyaGamma

////////////////////////////////////////////////////////////////////////////////
// Utility //
////////////////////////////////////////////////////////////////////////////////

void PolyaGamma::set_trunc(int trunc)
{

  if (trunc < 1) {
#ifndef NTHROW
    throw std::invalid_argument("PolyaGamma(int trunc): trunc < 1.");
#else
#ifndef USE_R
    fprintf(stderr, "Invalid parameter: PolyaGamma(int trunc): trunc < 1.  Setting trunc=1.\n");
#else
    Rprintf("Invalid parameter: PolyaGamma(int trunc): trunc < 1.  Setting trunc=1.\n");
#endif
    trunc = 1;
#endif
  }

  T = trunc;
  bvec.resize(T);

  for(int k=0; k < T; ++k){
    // + since we start indexing at 0.
    double d = ((double) k + 0.5);
    bvec[k] = FOURPISQ * d * d;
  }
} // set_trunc

double PolyaGamma::a(int n, double x)
{
  double K = (n + 0.5) * __PI;
  double y = 0;
  if (x > __TRUNC) {
    y = K * exp( -0.5 * K*K * x );
  }
  else if (x > 0) {
    double expnt = -1.5 * (log(0.5 * __PI)  + log(x)) + log(K) - 2.0 * (n+0.5)*(n+0.5) / x;
    y = exp(expnt);
  }
  return y;
}

double PolyaGamma::pigauss(double x, double Z)
{
  double b = sqrt(1.0 / x) * (x * Z - 1);
  double a = sqrt(1.0 / x) * (x * Z + 1) * -1.0;
  double y = R::pnorm(b,0,1,true, false) + exp(2 * Z) * R::pnorm(a,0,1,true, false);
  return y;
}

double PolyaGamma::mass_texpon(double Z, double fz)
{
  double t = __TRUNC;

  double b = sqrt(1.0 / t) * (t * Z - 1);
  double a = sqrt(1.0 / t) * (t * Z + 1) * -1.0;

  double x0 = log(fz) + fz * t;
  double xb = x0 - Z + R::pnorm(b,0,1,true, true);
  double xa = x0 + Z + R::pnorm(a,0,1,true, true);

  double qdivp = 4 / __PI * ( exp(xb) + exp(xa) );

  return 1.0 / (1.0 + qdivp);
}

double PolyaGamma::rtigauss(double Z)
{
  Z = fabs(Z);
  double t = __TRUNC;
  double X = t + 1.0;
  if (__TRUNC_RECIP > Z) { // mu > t
    double alpha = 0.0;
    while (R::runif(0,1) > alpha) {
      // X = t + 1.0;
      // while (X > t)
      // 	X = 1.0 / gamma_rate(0.5, 0.5);
      // Slightly faster to use truncated normal.
      double E1 = R::rexp(1.0); double E2 = R::rexp(1.0);
      while ( E1*E1 > 2 * E2 / t) {
        E1 = R::rexp(1.0); E2 = R::rexp(1.0);
      }
      X = 1 + E1 * t;
      X = t / (X * X);
      alpha = exp(-0.5 * X * (Z*Z));
    }
  }
  else {
    double mu = 1.0 / Z;
    while (X > t) {
      double Y = R::rnorm(0., 1.); Y *= Y;
      double half_mu = 0.5 * mu;
      double mu_Y    = mu  * Y;
      X = mu + half_mu * mu_Y - half_mu * sqrt(4 * mu_Y + mu_Y * mu_Y);
      if (R::runif(0,1) > mu / (mu + X))
        X = mu*mu / X;
    }
  }
  return X;
}

////////////////////////////////////////////////////////////////////////////////
// Sample //
////////////////////////////////////////////////////////////////////////////////

// double PolyaGamma::draw(double n, double z, RNG& r)
// {
//   return draw_sum_of_gammas(n, z, r);
// }

double PolyaGamma::draw(int n, double z)
{
  double sum = 0.0;
  for (int i = 0; i < n; ++i)
    sum += draw_like_devroye(z);
  return sum;
} // draw

double PolyaGamma::draw_sum_of_gammas(double n, double z)
{
  double x = 0;
  double kappa = z * z;
  for(int k=0; k < T; ++k)
    x += R::rgamma(n, 1.0) / (bvec[k] + kappa);
  return 2.0 * x;
} // draw_sum_of_gammas

double PolyaGamma::draw_like_devroye(double Z)
{
  // Change the parameter.
  Z = fabs(Z) * 0.5;

  // Now sample 0.25 * J^*(1, Z := Z/2).
  double fz = 0.125 * __PI*__PI + 0.5 * Z*Z;
  // ... Problems with large Z?  Try using q_over_p.
  // double p  = 0.5 * __PI * exp(-1.0 * fz * __TRUNC) / fz;
  // double q  = 2 * exp(-1.0 * Z) * pigauss(__TRUNC, Z);

  double X = 0.0;
  double S = 1.0;
  double Y = 0.0;
  // int iter = 0; If you want to keep track of iterations.

  while (true) {

    // if (unif() < p/(p+q))
    if ( R::runif(0,1) < mass_texpon(Z, fz) )
      X = __TRUNC + R::rexp(1.) / fz;
    else
      X = rtigauss(Z);

    S = a(0, X);
    Y = R::runif(0,1) * S;
    int n = 0;
    bool go = true;

    // Cap the number of iterations?
    while (go) {

      // Break infinite loop.  Put first so it always checks n==0.
      if (n % 1000 == 0) R_CheckUserInterrupt();

      ++n;
      if (n%2==1) {
        S = S - a(n, X);
        if ( Y<=S ) return 0.25 * X;
      }
      else {
        S = S + a(n, X);
        if ( Y>S ) go = false;
      }

    }
    // Need Y <= S in event that Y = S, e.g. when X = 0.

  }
} // draw_like_devroye

////////////////////////////////////////////////////////////////////////////////
// Static Members //
////////////////////////////////////////////////////////////////////////////////

double PolyaGamma::jj_m1(double b, double z)
{
  z = fabs(z);
  double m1 = 0.0;
  if (z > 1e-12)
    m1 = b * tanh(z) / z;
  else
    m1 = b * (1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6));
  return m1;
}

double PolyaGamma::jj_m2(double b, double z)
{
  z = fabs(z);
  double m2 = 0.0;
  if (z > 1e-12)
    m2 = (b+1) * b * pow(tanh(z)/z,2) + b * ((tanh(z)-z)/pow(z,3));
  else
    m2 = (b+1) * b * pow(1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6), 2) +
      b * ((-1.0/3) + (2.0/15) * pow(z,2) - (17.0/315) * pow(z,4));
  return m2;
}

double PolyaGamma::pg_m1(double b, double z)
{
  return jj_m1(b, 0.5 * z) * 0.25;
}

double PolyaGamma::pg_m2(double b, double z)
{
  return jj_m2(b, 0.5 * z) * 0.0625;
}
