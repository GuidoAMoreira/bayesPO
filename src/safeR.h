#ifndef __SAFER_H__
#define __SAFER_H__

/*
 * Code here is used to sample from the truncated normal distribution while being
 * thread safe.
 */

#ifndef M_SQRT_32
#define M_SQRT_32 5.656854249492380195206754896838792314276 // 128 bits to be safeR
#endif
#define M_1_SQRT_2PI 0.398942280401432677939946059934 // from R source
// #define M_SQRT2		1.414213562373095048801688724210 // from R source
#define M_2PI		6.283185307179586476925286766559 // from R source

#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

#ifndef ML_POSINF
#define ML_POSINF ((1.0) / 0.0)
#endif
#ifndef ML_NEGINF
#define ML_NEGINF ((-1.0) / 0.0)
#endif

#ifndef R_Q_P01_boundaries
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
if (log_p) {					                               \
  if(p == 0) /* upper bound*/			                \
return lower_tail ? _RIGHT_ : _LEFT_;	          \
  if(p == ML_NEGINF)				                        \
    return lower_tail ? _LEFT_ : _RIGHT_;	      \
}							                                        \
else { /* !log_p */					                        \
if(p == 0)					                                 \
  return lower_tail ? _LEFT_ : _RIGHT_;	        \
if(p == 1)					                                 \
  return lower_tail ? _RIGHT_ : _LEFT_;	        \
}
#endif

#ifndef R_D__0
#define R_D__0	(log_p ? ((-1.0) / 0.0) : 0.)		/* 0 */
#endif
#ifndef R_D__1
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#endif

#ifndef R_D_Lval
#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))
#endif
#ifndef R_DT_qIv
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) : R_D_Lval(p))
#endif

#ifndef R_D_Cval
#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))
#endif
#ifndef R_DT_CIv
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) : R_D_Cval(p))
#endif


// Shamelessly copied from R source.
// All thread unsafe code was (safely) removed.
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

static inline double safe_pnorm(double x, double mu, double sigma, double lower, double log_p) {
  double p, cp;
  p = (x - mu) / sigma;
  x = p;
  pnorm_both(x, &p, &cp, (lower ? 0 : 1), log_p);

  return(lower ? p : cp);
}

double safe_qnorm(double p, double mu, double sigma, int lower_tail, int log_p);

#endif
