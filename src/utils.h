#include "covariates.h"
#include "markovchain.h"
#include <RcppEigen.h>

// Create the mcStep object to be updated. Requires the pointer to a retrievCovs
// variable
mcStep startup_mcmc(Eigen::VectorXd beta, Eigen::VectorXd delta, double lambda,
                    std::string b_updater, std::string d_updater,
                    std::string l_updater, Rcpp::List parB, Rcpp::List parD,
                    Rcpp::List parL, std::string xClass,
                    Eigen::MatrixXd xValues, Eigen::VectorXi xI,
                    Eigen::VectorXi xO, retrievCovs *covars, double aD);

// Formating x functions
void importX_double(Eigen::MatrixXd x, int nb, int nd,
                    Eigen::VectorXi xI, Eigen::VectorXi xO,
                    std::vector<int> &x_data, Eigen::MatrixXd &zx_data,
                    Eigen::MatrixXd &wx_data);
void importX_int(Eigen::MatrixXi x, int nb, int nd,
                 Eigen::VectorXi xI, Eigen::VectorXi xO,
                 std::vector<int> &x_data, Eigen::MatrixXd &zx_data,
                 Eigen::MatrixXd &wx_data);
void determineX(Eigen::VectorXd x, int nb, int nd,
                      retrievCovs *b,
                      std::vector<int> &x_data, Eigen::MatrixXd &zx_data,
                      Eigen::MatrixXd &wx_data);

// Sampling. Simple inline function
inline Eigen::VectorXd sampleNormal(const Eigen::MatrixXd& Sigma)
{return Sigma.llt().matrixL() *
  Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rcpp::rnorm(Sigma.rows(), 0, 1));}
