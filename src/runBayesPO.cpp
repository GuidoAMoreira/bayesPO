#include "utils.h"
#include "covariates.h"
#include "markovchain.h"
#include "binary_regression.h"
#include "prior.h"
#include "link_functions.h"
#include <chrono>
#include <progress.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

enum LINKFUN {LOGIT, PROBIT};

// [[Rcpp::export]]
Rcpp::List runBayesPO(Eigen::VectorXd beta, Eigen::VectorXd delta,
                      double lambda,
                      Rcpp::List parB, Rcpp::List parD,
                      Rcpp::List parL, SEXP covariates,
                      int intensityLink, int observabilityLink, // Should be LINKFUN type, but it's giving error.
                      double areaD,
                      Eigen::MatrixXd xValues, // Observations covariates
                      Eigen::VectorXi intensityCovs, // Indices for intensity in background. Produced in R level.
                      Eigen::VectorXi observabilityCovs, // Indices for observability in background. Produced in R level.
                      Eigen::VectorXi xIntensityCovs, // Indices for intensity in observations. Produced in R level.
                      Eigen::VectorXi xObservabilityCovs, // Indices for observability in observations. Produced in R level.
                      int burnin, int thin, int iter, int threads, bool verbose)
{

  // Code is not thread safe yet.
  threads = 1;

#ifdef _OPENMP
  omp_set_num_threads( threads );
#endif

  Eigen::MatrixXd outBetas(iter/thin, beta.size()),
    outDeltas(iter/thin, delta.size());
  Eigen::VectorXd outLambdas(iter/thin), outLogPost(iter/thin),
    out_nU(iter/thin), out_nXp(iter/thin);
  std::vector<int> iC = std::vector<int>(
    intensityCovs.data(), intensityCovs.data() +
    intensityCovs.size()),
    oC = std::vector<int>(
      observabilityCovs.data(), observabilityCovs.data() +
        observabilityCovs.size());
  int i;


  // Let user know it's working
  if (verbose) Rcpp::Rcout << "Building MCMC model.\n";

  // Defining covariates object
  retrievCovs *covs = new retrievCovs_doubleMatrix(covariates, iC, oC);

  /*
   * double a, std::vector<int> x, Eigen::MatrixXd& zx,
   Eigen::MatrixXd& wx, retrievCovs* bg, double lambda,
   BinaryRegression* Int, BinaryRegression* Obs) : area(a),
   X(x), zX(zx), wX(wx), background(bg), lambdaStar(lambda),
   intensityRegression(Int), observabilityRegression(Obs)
   */

  // Create lean X matrices for efficient use.
  std::vector<int> x;
  Eigen::MatrixXd zx, wx;
  importX(xValues, intensityCovs.size(), observabilityCovs.size(),
          xIntensityCovs, xObservabilityCovs, x, zx, wx);

  // Get priors
  Eigen::VectorXd betaMean = parB["mean"];
  Eigen::MatrixXd betaCov = parB["covariance"];
  Eigen::VectorXd deltaMean = parD["mean"];
  Eigen::MatrixXd deltaCov = parD["covariance"];
  double lambdaA = parL["a"];
  double lambdaB = parL["b"];

  // Starting up MCMC
  LinkFunction *intLF, *obsLF;
  if (intensityLink) intLF = new Probit(zx); else intLF = new Logit(zx);
  if (observabilityLink) obsLF = new Probit(wx); else obsLF = new Logit(wx);
  MarkovChain mc(areaD, x, zx, wx, covs, lambdaA, lambdaB, lambda,
                 new LinearRegression(
                   new NormalPrior(betaMean, betaCov),
                   intLF, beta
                 ),
                 new LinearRegression(
                     new NormalPrior(deltaMean, deltaCov),
                     obsLF, delta
                 ));

  // Let user know it's working
  if (verbose) Rcpp::Rcout << "MCMC model built.\n";

  auto t1 = std::chrono::high_resolution_clock::now(); // Timing variable
  auto t2 = std::chrono::high_resolution_clock::now(); // Timing variable
  // Warming up the Markov Chain
  if (burnin)
  {
    if (verbose) Rcpp::Rcout << "Warming up the Markov Chain.\n";
    Progress progr_Burnin(burnin, true);
    t1 = std::chrono::high_resolution_clock::now();
    for (i = 0; i < burnin; i++)
      {
      progr_Burnin.increment();
      mc.update();
      }
    // Reset unObservedCounts after burn-in
    covs->resetUnobservedCounts();
    t2 = std::chrono::high_resolution_clock::now();
    if (verbose) Rcpp::Rcout << "Warm up complete. ";
    if (verbose) Rcpp::Rcout << "Warmup took " <<
      std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() <<
        " seconds.\n";
  }

  // Let user know it's working
  if (verbose) Rcpp::Rcout << "Sampling MCMC.\n";

#ifdef _OPENMP
  if ( threads > 0 )
    omp_set_num_threads( threads );
#endif

  Progress progr_Main(iter / thin, true);

  t1 = std::chrono::high_resolution_clock::now();
  // Actual MCMC
  for (i = 0; i < (iter/thin); i++)
  {
    mc.update(thin);
    /*    if (i == checkPointsIter[j])
     Rcpp::Rcout << "MCMC has completed " << 100*checkPoints[j++] <<
     "% of the iterations.\n";*/

    //R_CheckUserInterrupt();
    outBetas.row(i) = mc.getBeta();
    outDeltas.row(i) = mc.getDelta();
    outLambdas[i] = mc.getLambda();
    outLogPost[i] = mc.getLogPosterior();
    out_nU[i] = mc.getUsize();
    out_nXp[i] = mc.getXprimeSize();
    progr_Main.increment();
  }
  t2 = std::chrono::high_resolution_clock::now();
  //    auto t1 = std::chrono::high_resolution_clock::now();
  //    auto t2 = std::chrono::high_resolution_clock::now();
  //    Rcout << "Candidate choice took " << std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() << " microseconds.\n";

  // Let user know it's working
  if (verbose) Rcpp::Rcout << "MCMC sampling complete.\n";
  if (verbose) Rcpp::Rcout << "MCMC took " <<
    std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() <<
      " seconds.\n";

  return Rcpp::List::create(Rcpp::Named("beta") = outBetas,
                            Rcpp::Named("delta") = outDeltas,
                            Rcpp::Named("lambda") = outLambdas,
                            Rcpp::Named("nU") = out_nU,
                            Rcpp::Named("nXp") = out_nXp,
                            Rcpp::Named("logPost") = outLogPost,
                            Rcpp::Named("xPrimePred") = mc.getHeatMap());
}
