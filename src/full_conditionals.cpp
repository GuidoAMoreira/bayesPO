#include "full_conditionals.h"
using namespace arma;

double l_gamma(double *l, R_xlen_t x, R_xlen_t xp, R_xlen_t u, List p, double aD)
{
  double alpha = p["a"], beta = p["b"], out, parA = alpha + x + xp + u, parB = beta + aD;
  *l = R::rgamma(parA, 1 / parB);

  // log posterior contribution
  out = -(*l) * (parB) + std::log(*l)*parA;

  return out;
}

double bd_logit_normal(vec *bd, mat *oc, mat *zc, List p, double link(vec b, vec cov, bool c))
{
  unsigned long i, n1 = (*oc).n_rows, n0 = (*zc).n_rows;
  NumericVector b = p["mean"];
  NumericMatrix B = p["covariance"];

  PolyaGamma pg(1);
  mat V(&B[0],(*bd).size(),(*bd).size(),true,false);
  vec temp((*bd).size()), med = V * vec(b), xb1 = join_horiz(ones(n1),(*oc)) *
    (*bd), xb0 = join_horiz(ones(n0),(*zc)) * (*bd);

  // Calculating X' Omega X + B and X' kappa + B b
  for (i=0;i<n1;i++) // From the data matrix X
  {
    temp = join_vert(ones(1),(*oc).row(i).t());
    V += pg.draw_like_devroye(xb1[i]) * temp * temp.t();
    med += temp*0.5;
  }
  for (i=0;i<n0;i++) // From the data matrix X
  {
    temp = join_vert(ones(1),(*zc).row(i).t());
    V += pg.draw_like_devroye(xb0[i]) * temp * temp.t();
    med -= temp*0.5;
  }
  V = inv(V);
  *bd = chol(V) * randn((*bd).size()) + V * med;

  // log posterior contribution
  double out = 0;
  for (i=0;i<n1;i++)
  {
    out += link(*bd,(*oc).row(i).t(),false);
  }
  for (i=0;i<n0;i++)
  {
    out += link(*bd,(*zc).row(i).t(),true);
  }
  return out;
}
