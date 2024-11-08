#include "safeR.h"
#include <float.h>

// Shamelessly copied from R source
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p) {
  /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
   */
  const static double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const static double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };
  const static double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const static double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };
  const static double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const static double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };
  
  double xden, xnum, temp, del, eps, xsq, y;
  int i, lower, upper;
  
// #ifdef IEEE_754
//   if(ISNAN(x)) { *cum = *ccum = x; return; }
// #endif
  
  /* Consider changing these : */
  eps = DBL_EPSILON * 0.5;
  
  /* i_tail in {0,1,2} =^= {lower, upper, both} */
  lower = i_tail != 1;
  upper = i_tail != 0;
  
  y = fabs(x);
  if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
  if (y > eps) {
    xsq = x * x;
    xnum = a[4] * xsq;
    xden = xsq;
    for (i = 0; i < 3; ++i) {
      xnum = (xnum + a[i]) * xsq;
      xden = (xden + b[i]) * xsq;
    }
  } else xnum = xden = 0.0;
  
  temp = x * (xnum + a[3]) / (xden + b[3]);
  if(lower)  *cum = 0.5 + temp;
  if(upper) *ccum = 0.5 - temp;
  if(log_p) {
    if(lower)  *cum = log(*cum);
    if(upper) *ccum = log(*ccum);
  }
  }
  else if (y <= M_SQRT_32) {
    
    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
    
    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);
    
#define d_2(_x_) ldexp(_x_, -1) // == (_x_ / 2 )  "perfectly"
    
#define do_del(X)						                                   \
    xsq = ldexp(trunc(ldexp(X, 4)), -4);			               \
    del = (X - xsq) * (X + xsq);				                      \
    if(log_p) {						                                     \
      *cum = (-xsq * d_2(xsq)) -d_2(del) + log(temp);	    \
      if((lower && x > 0.) || (upper && x <= 0.))		       \
        *ccum = log1p(-exp(-xsq * d_2(xsq)) *		           \
          exp(-d_2(del)) * temp);		                       \
    }							                                              \
    else {							                                         \
      *cum = exp(-xsq * d_2(xsq)) * exp(-d_2(del)) * temp;\
      *ccum = 1.0 - *cum;					                            \
    }
  
#define swap_tail						                                \
  if (x > 0.) {/* swap  ccum <--> cum */			            \
    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
  }

do_del(y);
  swap_tail;
  }
  
  /* else	  |x| > sqrt(32) = 5.657 :
   * the next two case differentiations were really for lower=T, log=F
   * Particularly	 *not*	for  log_p !
   
   * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
   *
   * Note that we do want symmetry(0), lower/upper -> hence use y
   */
  else if((log_p && y < 1e170) /* avoid underflow below */
  /*  ^^^^^ MM FIXME: could speed up for log_p and  y := |x| >> 5.657 !
   * Then, make use of  Abramowitz & Stegun, 26.2.13, p.932,  something like
   
   * Even smarter: work with   example(pnormAsymp, package="DPQ")
   
   xsq = x*x;
   
   if(xsq * DBL_EPSILON < 1.)
   del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
   else
   del = 0.;
   *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
   *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
   
   swap_tail;
   
   Yes, but xsq might be infinite;
   well, actually  x = -1.34..e154 = -sqrt(DBL_MAX) already overflows x^2
   The largest x for which  x/2*x is finite is
   x = +/- 1.89615038e154 ~= sqrt(2) * sqrt(.Machine$double.xmax)
   
   NB: allowing "DENORMS" ==> boundaries at +/- 38.4674  <--> qnorm(log(2^-1074), log.p=TRUE)
   --                               rather than 37.5193 (up to R 4.4.x)
   */
  || (lower && -38.4674 < x  &&  x < 8.2924)
            || (upper && -8.2924  < x  &&  x < 38.4674)
  ) {
    
    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
    xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (M_1_SQRT_2PI - temp) / y;
    
    do_del(x);
    swap_tail;
  } else { /* large |x| such that probs are 0 or 1 */
    if(x > 0) {	*cum = R_D__1; *ccum = R_D__0;	}
    else {	        *cum = R_D__0; *ccum = R_D__1;	}
  }

  return;
}

double qnorm(double p, double mu, double sigma, int lower_tail, int log_p)
{
  double p_, q, r, val;
  
  p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;
    
    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     
     Produces the normal deviate Z corresponding to a given lower
     tail area of P; Z is accurate to about 1 part in 10**16.
     
     (original fortran code used PARAMETER(..) for the coefficients
     and provided hash codes for checking them...)
     */
    if (fabs(q) <= .425) {/* |p~ - 0.5| <= .425  <==> 0.075 <= p~ <= 0.925 */
    r = .180625 - q * q; // = .425^2 - q^2  >= 0
      val =
        q * (((((((r * 2509.0809287301226727 +
        33430.575583588128105) * r + 67265.770927008700853) * r +
        45921.953931549871457) * r + 13731.693765509461125) * r +
        1971.5909503065514427) * r + 133.14166789178437745) * r +
        3.387132872796366608)
        / (((((((r * 5226.495278852854561 +
          28729.085735721942674) * r + 39307.89580009271061) * r +
          21213.794301586595867) * r + 5394.1960214247511077) * r +
          687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary :
     *  r := log(p~);  p~ = min(p, 1-p) < 0.075 :  */
    double lp;
      if(log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) {
        lp = p;
      } else {
        lp = log( (q > 0) ? R_DT_CIv(p) /* 1-p */ : p_ /* = R_DT_Iv(p) ^=  p */);
      }
      // r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 ) :
      r = sqrt(-lp);
      if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
    r += -1.6;
        val = (((((((r * 7.7454501427834140764e-4 +
                          .0227238449892691845833) * r + .24178072517745061177) *
        r + 1.27045825245236838258) * r +
        3.64784832476320460504) * r + 5.7694972214606914055) *
        r + 4.6303378461565452959) * r +
        1.42343711074968357734)
          / (((((((r *
            1.05075007164441684324e-9 + 5.475938084995344946e-4) *
            r + .0151986665636164571966) * r +
                .14810397642748007459) * r + .68976733498510000455) *
            r + 1.6763848301838038494) * r +
            2.05319162663775882187) * r + 1.);
      }
      else if(r <= 27) { /* p is very close to  0 or 1: r in (5, 27] :
     *  r >   5 <==> min(p,1-p)  < exp(-25) = 1.3888..e-11
     *  r <= 27 <==> min(p,1-p) >= exp(-27^2) = exp(-729) ~= 2.507972e-317
     * i.e., we are just barely in the range where min(p, 1-p) has not yet underflowed to zero.
     */
    // Wichura, p.478: minimax rational approx R_3(t) is for 5 <= t <= 27  (t :== r)
    r += -5.;
        val = (((((((r * 2.01033439929228813265e-7 +
          2.71155556874348757815e-5) * r +
           .0012426609473880784386) * r + .026532189526576123093) *
          r + .29656057182850489123) * r +
          1.7848265399172913358) * r + 5.4637849111641143699) *
          r + 6.6579046435011037772)
          / (((((((r *
            2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
            r + 1.8463183175100546818e-5) * r +
            7.868691311456132591e-4) * r + .0148753612908506148525)
                 * r + .13692988092273580531) * r +
                       .59983220655588793769) * r + 1.);
      }
      else { // r > 27: p is *really* close to 0 or 1 .. practically only when log_p =TRUE
        if(r >= 6.4e8) { // p is *very extremely* close to 0 or 1
          // Using the asymptotical formula ("0-th order"): qn = sqrt(2*s)
          val = r * M_SQRT2;
        } else {
          double s2 = -ldexp(lp, 1), // = -2*lp = 2s
            x2 = s2 - log(M_2PI * s2); // = xs_1
          // if(r >= 36000.)  # <==> s >= 36000^2   use x2 = xs_1  above
          if(r < 36000.) {
            x2 = s2 - log(M_2PI * x2) - 2./(2. + x2); // == xs_2
            if(r < 840.) { // 27 < r < 840
              x2 = s2 - log(M_2PI * x2) + 2*log1p(- (1 - 1/(4 + x2))/(2. + x2)); // == xs_3
              if(r < 109.) { // 27 < r < 109
                x2 = s2 - log(M_2PI * x2) +
                  2*log1p(- (1 - (1 - 5/(6 + x2))/(4. + x2))/(2. + x2)); // == xs_4
                if(r < 55.) { // 27 < r < 55
                  x2 = s2 - log(M_2PI * x2) +
                    2*log1p(- (1 - (1 - (5 - 9/(8. + x2))/(6. + x2))/(4. + x2))/(2. + x2)); // == xs_5
                }
              }
            }
          }
          val = sqrt(x2);
        }
      }
      if(q < 0.0)
        val = -val;
    }
    return mu + sigma * val;
}

