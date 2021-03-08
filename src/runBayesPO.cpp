#include "auxiliary_functions.h"
#include "covariates.h"
#include "markovchain.h"
#include <chrono>
using namespace Rcpp;

// [[Rcpp::export]]
List runBayesPO(NumericVector beta, NumericVector delta, double lambda,
                String b_updater, String d_updater, String l_updater, List parB,
                List parD, List parL, String covsClass, SEXP covariates,
                double areaD, String xClass, SEXP xValues,
                IntegerVector intensityCovs, IntegerVector observabilityCovs,
                IntegerVector xIntensityCovs, IntegerVector xObservabilityCovs,
                int burnin, int thin, int iter)
{
  NumericMatrix outBetas(iter/thin,beta.size()), outDeltas(iter/thin,delta.size());
  NumericVector outLambdas(iter/thin), outLogPost(iter/thin), out_nU(iter/thin), out_nXp(iter/thin);
  R_xlen_t i;
  arma::vec checkPoints({0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9});
  arma::vec::iterator checked = checkPoints.begin();
  R_xlen_t stopIter = R_xlen_t(checked[0] * iter);
  auto t1 = std::chrono::high_resolution_clock::now(); // Timing variable
  auto t2 = std::chrono::high_resolution_clock::now(); // Timing variable

  // Let user know it's working
  Rcout << "Building MCMC model.\n";

  // Defining covariates object
  retrievCovs *covs;
  if (covsClass == "num_mat")
  {
    covs = new retrievCovs_doubleMatrix(covariates,
                                        as<std::vector<R_xlen_t> >(intensityCovs),
                                        as<std::vector<R_xlen_t> >(observabilityCovs));
  } else if (covsClass == "int_mat")
  {
    covs = new retrievCovs_intMatrix(covariates,
                                     as<std::vector<R_xlen_t> >(intensityCovs),
                                     as<std::vector<R_xlen_t> >(observabilityCovs));
  } else if (covsClass == "std_normal") {covs = new retrievCovs_normal(
    as<std::vector<R_xlen_t> >(intensityCovs),
    as<std::vector<R_xlen_t> >(observabilityCovs),beta.size()-1,delta.size()-1);}

  // Starting up MCMC
  mcStep MarkovChain = startup_mcmc(beta,delta,lambda,b_updater,d_updater,
                                    l_updater,parB,parD,parL,xClass,xValues,
                                    xIntensityCovs,xObservabilityCovs,covs,areaD);

  // Let user know it's working
  Rcout << "MCMC model built.\n";

  // Warming up the Markov Chain
  if (burnin)
    {
    Rcout << "Warming up the Markov Chain.\n";
    stopIter = R_xlen_t(checked[0] * burnin);
    t1 = std::chrono::high_resolution_clock::now();
    for (i=0;i<burnin;i++)
      {
      MarkovChain.update();
      if (i == stopIter)
      {
        Rcout << "Warmup has completed " << 100*checked[0] << "% of the iterations.\n";
        stopIter = R_xlen_t((++checked)[0] * burnin);
      }
      }
    t2 = std::chrono::high_resolution_clock::now();
    Rcout << "Warm up complete. ";
    Rcout << "Warmup took " <<
      std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << " seconds.\n";
    checked = checkPoints.begin();
    stopIter = R_xlen_t(checked[0] * iter);
    }

  // Let user know it's working
  Rcout << "Sampling MCMC.\n";

  t1 = std::chrono::high_resolution_clock::now();
  // Actual MCMC
  for (i=0;i<iter;i++)
  {
    MarkovChain.update();
    if (i == stopIter)
    {
      Rcout << "MCMC has completed " << 100*checked[0] << "% of the iterations.\n";
      stopIter = R_xlen_t((++checked)[0] * iter);
    }
    if (i % 1000) {R_CheckUserInterrupt();}
    if (i % thin == 0) // Storing results
    {
      outBetas(i/thin,_) = NumericVector(MarkovChain.beta.begin(),MarkovChain.beta.end());
      outDeltas(i/thin,_) = NumericVector(MarkovChain.delta.begin(),MarkovChain.delta.end());
      outLambdas[i/thin] = MarkovChain.lambdaStar;
      outLogPost[i/thin] = MarkovChain.logPosterior;
      out_nU[i/thin] = MarkovChain.U.size();
      out_nXp[i/thin] = MarkovChain.Xprime.size();
    }
  }
  t2 = std::chrono::high_resolution_clock::now();
  //    auto t1 = std::chrono::high_resolution_clock::now();
  //    auto t2 = std::chrono::high_resolution_clock::now();
  //    Rcout << "Candidate choice took " << std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() << " microseconds.\n";

  // Let user know it's working
  Rcout << "MCMC sampling complete.\n";
  Rcout << "MCMC took " <<
    std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << " seconds.\n";

  List output = List::create(Named("beta") = outBetas, Named("delta") = outDeltas,
                             Named("lambda") = outLambdas, Named("nU") = out_nU,
                             Named("nXp") = out_nXp, Named("logPost") = outLogPost);
  return output;
}

