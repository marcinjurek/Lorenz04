#include "vectorTools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;




double Wn(const arma::colvec& X, const int& n, const int& K){

  int N = X.n_rows;
  int J = K / 2;

  double first = X(mod(n - J, N)) / K;
  double last = X(mod(n + J, N)) / K;
  if( K % 2 == 0){
    first /= 2;
    last  /= 2;
  }
  
  double sums = first + last - X(mod(n, N)) / K;
  
  for(int i = 0; i < J; ++i)  {
    sums = sums + X(mod(n - i, N)) / K + X(mod(n + i, N)) / K;
  }

  return(sums);
  
}




double XX(const arma::colvec& X, const int& n, const int& K){
  
  int N = X.n_rows;
  int J = K / 2;
  
  double first = (Wn(X, mod(n - K - J, N), K) * X(mod(n + K - J, N))) / K;
  double last  = (Wn(X, mod(n - K + J, N), K) * X(mod(n + K + J, N))) / K;

  if( K % 2 == 0){
    first /= 2;
    last /= 2;
  }

  double r_sum = first + last + (Wn(X, mod(n - K, N), K) * X(mod(n + K, N)) / K);

  for(int i = 1; i < J; ++i){
    double jcontrib = Wn(X, mod(n - K + i, N), K) * X(mod(n + K + i, N)) / K;
    double negjcontrib = Wn(X, mod(n - K - i, N), K) * X(mod(n + K - i, N)) / K;
    r_sum = r_sum + jcontrib + negjcontrib;
  }
  
  double XX_val = (-Wn(X, mod(n - (2 * K), N), K)) * Wn(X, mod(n - K, N), K) + r_sum;

  if(XX_val==arma::datum::nan){
    Rcout << "ERROR!" << std::endl;
  }
  
  return XX_val;
}





arma::vec getXFromZ(const arma::vec& Z, const double& alpha, const double& beta, const int& I){

  int N = Z.n_rows;
  arma::vec X = arma::zeros<arma::vec>(N);
  for(int j = 0; j < N; ++j){

    double first  = Z(mod(j + I, N)) * (alpha - beta*I)/2;
    double last   = Z(mod(j - I, N)) * (alpha - beta*I)/2;
    double middle = Z(j) * alpha;

    double Xj = first + last + middle;

    for(int i = 1; i < I; ++i){
      Xj += (Z(mod(j + i, N)) + Z(mod(j - i, N))) * (alpha - beta*i);
    }
    X(j) = Xj;
  }

  return X;
}





double bracket(const arma::colvec& X, const arma::colvec& Y, const int& K, const int& n){

  int N = X.n_rows;
  
  //if X==Y use a faster algorithm from eq. 10  
  if( K==1 ){
    if(n==0){
      //Rcout << "here" << std::endl;
    }
    return -X(mod(n - 2, N)) * Y(mod(n - 1, N)) + X(mod(n - 1, N)) * Y(mod(n + 1, N));
  } else if(sum(abs(X - Y)) < 1e-10){
    if(n==0){
      //Rcout << "not here" << std::endl;
    }
    return XX(X, n, K);
  } else {
    stop("Inadmissible value of K");
    return -1;
  }   
  
}




arma::vec scalarRHS(const arma::vec& Z, const int& K, const double& F, const int& I,
		    const double& b, const double& c, const double& alpha, const double& beta){

  int N = Z.n_rows;
  
  if( Z.has_nan()){
    stop("nan's detected in Z");
  }
  
  arma::vec X = getXFromZ(Z, alpha, beta, I);
  arma::vec Y = Z - X;

  if(X.has_nan() | Y.has_nan()){
    stop("nan's detected in X or Y");
  }

  arma::vec k = arma::zeros<arma::vec>(N);
  
  for(int j = 0; j < N; ++j){
    double bXX = bracket(X, X, K, j);
    double bYY = bracket(Y, Y, 1, j);
    double bYX = bracket(Y, X, 1, j);
   
    k(j) = bXX + b*b*bYY + c*bYX - X(j) - b*Y(j) + F;
    
  }
  
  if(k.has_nan()){
    stop("nan's detected in k");
  }
  
  return k;

}
