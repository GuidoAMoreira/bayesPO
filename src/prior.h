#ifndef __PRIOR_BAYESPO_H__
#define __PRIOR_BAYESPO_H__

#include <RcppEigen.h>

#ifndef M_LOG_SQRT2PI
#define M_LOG_SQRT2PI 0.9189385332046727222894705738087400271472
#endif

/*
 * WARNING: This code has not been thread safe proofed. Do NOT attempt to
 * run with multiple cores.
 */

class Prior{
public:
  virtual double logDensity(const Eigen::VectorXd&) = 0;
  virtual Eigen::VectorXd sampleMP(const Eigen::VectorXd& meanVec,
                                   const Eigen::MatrixXd& covMat) = 0;
  virtual ~Prior() = default;
};

class NormalPrior : public Prior {
  const Eigen::VectorXd m;
  const Eigen::MatrixXd v;
  const Eigen::LLT<Eigen::MatrixXd> vChol;
  const Eigen::MatrixXd p;
  const Eigen::VectorXd pm;
  const double halflog_pDet;

  Eigen::LLT<Eigen::MatrixXd> solver;
public:
  double logDensity(const Eigen::VectorXd& vec) {
    Eigen::VectorXd xMinusm = vec - m;
    return -M_LOG_SQRT2PI + halflog_pDet -
      0.5 * xMinusm.transpose() * p * xMinusm;
  }
  Eigen::VectorXd sampleMP(const Eigen::VectorXd& meanVec,
                           const Eigen::MatrixXd& covMat) {
    solver.compute(covMat + p);
    Eigen::VectorXd stdNormal;

#pragma omp critical
    stdNormal = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rcpp::rnorm(v.rows(), 0, 1));

    return solver.matrixU().transpose().solve(stdNormal) +
      solver.solve(meanVec + pm);
  }

  NormalPrior(Eigen::VectorXd mu, Eigen::MatrixXd Sigma) :
    m(mu), v(Sigma), vChol(Sigma.llt()),
    p(vChol.solve(Eigen::MatrixXd::Identity(Sigma.rows(), Sigma.cols()))),
    pm(p * m), halflog_pDet(vChol.matrixLLT().diagonal().array().log().sum()) {}
};


#endif
