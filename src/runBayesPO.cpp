#include "utils.h"
#include "covariates.h"
#include "markovchain.h"
#include <chrono>
#include <progress.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::List runBayesPO(Eigen::VectorXd beta, Eigen::VectorXd delta,
                      double lambda, Rcpp::String b_updater,
                      Rcpp::String d_updater, Rcpp::String l_updater,
                      Rcpp::List parB, Rcpp::List parD,
                      Rcpp::List parL, Rcpp::String covsClass, SEXP covariates,
                      double areaD, Rcpp::String xClass,
                      Eigen::MatrixXd xValues,
                      Eigen::VectorXi intensityCovs,
                      Eigen::VectorXi observabilityCovs,
                      Eigen::VectorXi xIntensityCovs,
                      Eigen::VectorXi xObservabilityCovs,
                      int burnin, int thin, int iter, int threads, bool verbose)
{

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


  auto t1 = std::chrono::high_resolution_clock::now(); // Timing variable
  auto t2 = std::chrono::high_resolution_clock::now(); // Timing variable

  // Let user know it's working
  if (verbose) Rcpp::Rcout << "Building MCMC model.\n";

  // Defining covariates object
  retrievCovs *covs;
  if (covsClass == "num_mat")
  {
    covs = new retrievCovs_doubleMatrix(covariates, iC, oC);
  } else if (covsClass == "int_mat")
  {
    covs = new retrievCovs_intMatrix(covariates, iC, oC);
  } else if (covsClass == "std_normal") {covs = new retrievCovs_normal(iC, oC);}
  else {Rcpp::stop("Unidentified background data type - C++ level");}

  // Starting up MCMC
  mcStep MarkovChain = startup_mcmc(beta, delta, lambda, b_updater, d_updater,
                                    l_updater, parB, parD, parL, xClass, xValues,
                                    xIntensityCovs, xObservabilityCovs, covs,
                                    areaD);

  // Let user know it's working
  if (verbose) Rcpp::Rcout << "MCMC model built.\n";

  // Warming up the Markov Chain
  if (burnin)
  {
    if (verbose) Rcpp::Rcout << "Warming up the Markov Chain.\n";
    Progress progr_Burnin(burnin, true);
    t1 = std::chrono::high_resolution_clock::now();
    for (i = 0; i < burnin; i++)
      {
      progr_Burnin.increment();
      MarkovChain.update();
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
    MarkovChain.update(thin);
    /*    if (i == checkPointsIter[j])
     Rcpp::Rcout << "MCMC has completed " << 100*checkPoints[j++] <<
     "% of the iterations.\n";*/

    //R_CheckUserInterrupt();
    outBetas.row(i) = MarkovChain.beta->effects;
    outDeltas.row(i) = MarkovChain.delta->effects;
    outLambdas[i] = MarkovChain.lambda->l;
    outLogPost[i] = MarkovChain.logPosterior;
    out_nU[i] = MarkovChain.U.size();
    out_nXp[i] = MarkovChain.Xprime.size();
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
                            Rcpp::Named("xPrimePred") =
                              MarkovChain.background->getUnobservedCounts());
}
