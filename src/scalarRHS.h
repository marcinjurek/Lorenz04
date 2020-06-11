#include "vectorTools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec scalarRHS(const arma::vec& Z, const int& K, const double& F, const int& I, const double& b, const double& c, const double& alpha, const double& beta);
