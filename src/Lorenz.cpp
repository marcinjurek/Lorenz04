#include "scalarRHS.h"
#include "vectorRHS.h"
#include <RcppArmadillo.h>
#include "exactGrad.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;





//' A single iteration of the Lorenz 04 model
//'
//' @param Z0 initial state
//' @param F F from the Lorenz 04 model
//' @param K K from the Lorenz 04 model
//' @param dt time step
//' @param M the number of steps that make a time step
//' @param vectorAlgo whether to run the vector version of the algorithm
//' @return dx the increment to reach the new state
//' @export
// [[Rcpp::export]]
arma::vec DeltaLorenz04M2Cpp(const arma::vec& Z0, const double & F, const int& K,
                             const double& dt, const int& M, const bool& vectorAlgo){


  int I = 12;
  double b = 10;
  double c = 2.5;
  //double c = 0.6;
  double alpha = (3*I*I + 3.0)/(2*I*I*I + 4*I);
  double beta = (2*I*I + 1.0)/(I*I*I*I + 2*I*I);

  
  // CapInitialize vector of zeros that will "update" XX
  arma::vec Z = Z0;

  arma::vec k1, k2, k3, k4;


  for(int i = 0; i < M ; ++i){

    if(vectorAlgo){
      
      k1 = vectorRHS(Z, K, F);
      k2 = vectorRHS(Z + 0.5*dt*k1, K, F);
      k3 = vectorRHS(Z + 0.5*dt*k2, K, F);
      k4 = vectorRHS(Z + dt*k3, K, F);
      
    } else {

      k1 = scalarRHS(Z, K, F, I, b, c, alpha, beta);
      k2 = scalarRHS(Z + 0.5*dt*k1, K, F, I, b, c, alpha, beta);
      k3 = scalarRHS(Z + 0.5*dt*k2, K, F, I, b, c, alpha, beta);
      k4 = scalarRHS(Z + dt*k3, K, F, I, b, c, alpha, beta);
      
    }

    arma::vec dz = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
    Z += dz;
      
  }
  
  return Z - Z0;
}


//' Simulate from the Lorenz 04 model
//'
//' @param Xinit initial state
//' @param F_Lor F from the Lorenz 04 model
//' @param K_Lor K from the Lorenz 04 model
//' @param dt time step
//' @param M the number of steps that make a time step
//' @param iter number of iterations
//' @param burn how many steps to discard
//' @param vectorAlgo whether to run the vector version of the algorithm
//' @return Xiter values of all iterations
//' @export
// [[Rcpp::export]]
arma::mat Lorenz04M2SimCpp(const arma::vec& Xinit, const int& F_Lor, const int& K_Lor,
                           const double& dt, const int& M, const int& iter, const int& burn,
			   const bool& vectorAlgo){
  // Get N_Lor from Xinit
  int N_Lor = Xinit.n_rows;
  // Create Matrix to store generated data
  arma::mat Xiter(N_Lor, iter+1);
  
  // Create vector to allow for burn-in iterations
  arma::vec newXburn = Xinit;
  
  int burnin = burn;
  
  if( burnin > 0 ){
    for(int i = 0; i < burnin; ++i){
      arma::vec delta = DeltaLorenz04M2Cpp(newXburn, F_Lor, K_Lor, dt, M, vectorAlgo);
      newXburn = newXburn + delta;
    }
  } 
  
  Xiter.col(0) = newXburn;
  
  for(int i = 1; i < iter+1; ++i){
    arma::vec delta = DeltaLorenz04M2Cpp(Xiter.col(i - 1), F_Lor, K_Lor, dt, M, vectorAlgo);
    Xiter.col(i) = Xiter.col(i-1) + delta;
  }
  
  return Xiter.tail_cols(iter);
}





arma::mat EGradient(const arma::vec& X, const int& K, const double& dt, const double& F){
  
  int N = X.n_rows;
  
  arma::vec Xp = X + 0.5*dt*vectorRHS(X, K, F);
  arma::vec Xb = X + 0.5*dt*vectorRHS(Xp, K, F);
  arma::vec Xt = X + dt*vectorRHS(Xb, K, F);
  
  arma::mat A = exactGrad(X, K);
  arma::mat B = exactGrad(Xp, K) % (1 + 0.5*dt*A);
  arma::mat C = exactGrad(Xb, K) % (1 + 0.5*dt*B);
  arma::mat D = exactGrad(Xt, K) % (1 + dt*C);
  
  return (A + B + C + D)*(dt/6) + arma::eye(N,N);
}



//' Gradient for the Lorenz model
//'
//' @param X initial state
//' @param K K from the Lorenz 04 model
//' @param M the number of steps that make a time step
//' @param dt time step
//' @param F F from the Lorenz 04 model
//' @export
// [[Rcpp::export]]
arma::mat exactGradient(const arma::vec& X, const int& K, const int& M, const double& dt, const double& F){

  int N = X.n_rows;
  
  if(K>=(N/2)){
    stop("K is too big for this N.");
  }
  
  arma::mat G = arma::eye(N, N);
  arma::vec Xnew = X;
  for( int m=1; m<=M; m++){
    G = EGradient(Xnew, K, dt, F) * G;
    Xnew = Xnew + DeltaLorenz04M2Cpp(Xnew, F, K, dt, 1, false);
  }
  return( G );
}




arma::mat Lorenz04M2Gradient(const arma::vec& X, const int& F_Lor, const int& K_Lor,
                           const double& dt, const int& M, const int& iter, const int& burn,
			     const bool& newAlgo){

  
  int N = X.n_rows;
  arma::vec newf(N);
  arma::vec del(N);
  arma::vec tempX = X;
  arma::vec reff = X + DeltaLorenz04M2Cpp(X, F_Lor, K_Lor, dt, M, newAlgo);
  arma::mat gradient = arma::zeros<arma::mat>(N,N);

  arma::vec pert = arma::ones<arma::vec>(N)*1e-8;
  arma::vec delta = arma::max(pert, arma::abs(X % pert));
  
  
  for(arma::uword i=0; i<N; i++){
  
    tempX(i) = tempX(i) + delta(i);
    newf = tempX + DeltaLorenz04M2Cpp(tempX, F_Lor, K_Lor, dt, M, newAlgo);
    del = (newf - reff) / delta(i);
    
    gradient.col(i) = del;
    tempX(i) = X(i);
  }
  
  return gradient;

}
