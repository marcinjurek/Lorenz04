#include <RcppArmadillo.h>
#include "vectorTools.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



double sumPrime(arma::vec v, const int& K, const int& p, const int& ind1, const int& ind2){
  if( (K % 2)==0 ){
    v.head(1) /=2;
    if(v.size()>1){
      v.tail(1) /= 2;
    }
  }
  return(sum(v));
}



double getMultiplier(const int& K, const int& i, const int& j, const int& J){
  
  double multiplier = 1;
  if( (K % 2)==0 ){
    if( abs(j)==J ){
      multiplier /= 2;
    }
    if( abs(i)==J ){
      multiplier /= 2;
    } 
  }
  return multiplier;
}



bool check(int ind1, int ind2, int p){
  if(ind1 < ind2){
    return p>=ind1 && p<=ind2;
  } else if(ind1==ind2) {
    return p==ind1;
  } else {
    return p>=ind1 || p<=ind2;
  }
}


bool isIn1( int n, int p, int K, int J, int N ){
  int ind1 = mod(n - 2*K - J, N);
  int ind2 = mod(n - 2*K + J, N);
  return check(ind1, ind2, p);
}


bool isBorder1( int n, int p, int K, int J, int N ){
  if( K%2==1 ){
    return 0;
  }
  int ind1 = mod(n - 2*K - J, N);
  int ind2 = mod(n - 2*K + J, N);
  return( p==ind1 || p == ind2 );
}


bool isIn2( int n, int p, int K, int J, int N ){
  int ind1 = mod(n - K - J, N);
  int ind2 = mod(n - K + J, N);
  return check(ind1, ind2, p);
}


bool isBorder2( int n, int p, int K, int J, int N ){
  if( K%2==1 ){
    return 0;
  }
  int ind1 = mod(n - K - J, N);
  int ind2 = mod(n - K + J, N);
  return( p==ind1 || p==ind2 );
}


bool isIn3( int n, int p, int K, int J, int N ){
  int ind1 = mod(n - K - 2*J, N);
  int ind2 = mod(n - K + 2*J, N);
  return check(ind1, ind2, p);
}


bool isIn4( int n, int p, int K, int J, int N ){
  int ind1 = mod(n + K - J, N);
  int ind2 = mod(n + K + J, N);
  return check(ind1, ind2, p);
}


bool isBorder4( int n, int p, int K, int J, int N ){
  if( K%2==1 ){
    return 0;
  }
  int ind1 = mod(n + K - J, N);
  int ind2 = mod(n + K + J, N);
  return( p == ind1 || p == ind2);
}



//' Calculate the exact gradient in the Lorenz model
//'
//' @param X the vector at which to evaluate the gradient
//' @param K parameter from the Lorenz model
//' @return grad the matrix with the gradient
//' @export
// [[Rcpp::export]]
arma::mat exactGrad(const arma::vec& X, const int& K) {
  
  int N = X.n_rows;
  mat grad(N,N);
  int J = K/2;
  
  uword ind1, ind2;
  
  for( int n = 0; n<N; n++ ){
    for( int p = 0; p<N; p++ ){
      double d = 0;
      if( isIn1(n, p, K, J, N) ){

        ind1 = mod(n - K - J, N);
        ind2 = mod(n - K + J, N);
        vec masked = mask(X, ind1, ind2, true);
        masked *= (1-isBorder1(n, p, K, J, N) * 0.5);
        d -= sumPrime(masked, K, p, ind1, ind2);
      }
      if( isIn2(n, p, K, J, N) ){

        ind1 = mod(n - 2*K - J, N);
        ind2 = mod(n - 2*K + J, N);
        vec masked = mask(X, ind1, ind2, true);
        masked *= (1-isBorder2(n, p, K, J, N) * 0.5);
        d -= sumPrime(masked, K, p, ind1, ind2); 
      }
      if( isIn3(n, p, K, J, N) ){
        for( int j=-J; j<=J; j++){
          for( int i=-J; i<=J; i++){
            if(mod(n - K + j - i, N) == p){
              d += X(mod(n + K + j, N)) * getMultiplier(K, i, j, J);
            }
          }
        }
      }
      if( isIn4(n, p, K, J, N) ){
        int j = p - n - K;
        ind1 = mod(n - K + j - J, N);
        ind2 = mod(n - K + j + J, N);
        vec masked = mask(X, ind1, ind2, true);
        masked *= (1-isBorder4(n, p, K, J, N) * 0.5);
        d += sumPrime(masked, K, p, ind1, ind2);
      }
      if( p==n ){
        d -= 1;
      }
      grad(n, p) = d/(K*K);
    }
  }
  return grad;
}

