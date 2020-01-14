#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


int mod_cpp(const int& i, const int& n){ 
  int mod = (i - 1) % n;
  if(mod < 0)
  { mod += n; }
  return mod;
}


double Wn_Even_Cpp(const arma::colvec& XX, const int& n, const int& k, const int& N){
  
  double K = k;
  double J = K / 2;
  
  double first = XX(mod_cpp(n - J, N)) / (2 * K);
  double last = XX(mod_cpp(n + J, N)) / (2 * K);
  
  double sums = first + last - XX(mod_cpp(n, N)) / K;
  
  for(int i = 0; i < J; ++i)
  {
    sums = sums + XX(mod_cpp(n - i, N)) / K + XX(mod_cpp(n + i, N)) / K;
  }
  
  return(sums);
}


double Wn_Odd_Cpp(const arma::colvec& XX, const int& n, const int& k, const int& N){
  double K = k;
  double J = (K - 1) /2;
  double sums = - XX(mod_cpp(n, N)) / K;
  
  for(int i = 0; i <= J; ++i)
  {
    sums = sums + XX(mod_cpp(n - i, N)) / K + XX(mod_cpp(n + i, N)) / K;
  }
  
  return(sums);
}


double XX_Kn_Even_Cpp(const arma::colvec& XX, const int& n, const int& k){
  
  int N = XX.n_rows;
  double K = k;
  double J = K / 2;
  
  double first = (Wn_Even_Cpp(XX, mod_cpp(n - K - J, N), K, N) * XX(mod_cpp(n + K - J, N))) / (2 * K);
  double last = (Wn_Even_Cpp(XX, mod_cpp(n - K + J, N), K, N) * XX(mod_cpp(n + K + J, N))) / (2 * K);
  
  double r_sum = first + last + (Wn_Even_Cpp(XX, mod_cpp(n - K, N), K, N) * XX(mod_cpp(n + K, N)) / K);
  
  for(int i = 1; i < J; ++i){
    r_sum = r_sum + Wn_Even_Cpp(XX, mod_cpp(n - K - i, N), K, N) * XX(mod_cpp(n + K - i, N)) / K + Wn_Even_Cpp(XX, mod_cpp(n - K + i, N), K, N) * XX(mod_cpp(n + K + i, N)) / K;
    //std::cout << "iteration " << i << ": " << r_sum << std::endl;
  }
  
  double XX_Kn_val = (-Wn_Even_Cpp(XX, mod_cpp(n - (2 * K), N), K, N)) * Wn_Even_Cpp(XX, mod_cpp(n - K, N), K, N) + r_sum;
  //std::cout << (-Wn_Even_Cpp(XX, mod_cpp(n - (2 * K), N), K, N)) * Wn_Even_Cpp(XX, mod_cpp(n - K, N), K, N) << std::endl;
  return XX_Kn_val;
}


double XX_Kn_Odd_Cpp(const arma::colvec& XX, const int& n, const int& k){
  
  int N = XX.n_rows;
  double K = k;
  double J = (K - 1) / 2;
  
  double r_sum = (Wn_Odd_Cpp(XX, mod_cpp(n - K, N), K, N) * XX(mod_cpp(n + K, N)) / K);
  
  for(int i = 1; i <= J; ++i){
    r_sum = r_sum + Wn_Odd_Cpp(XX, mod_cpp(n - K - i, N), K, N) * XX(mod_cpp(n + K - i, N)) / K + Wn_Odd_Cpp(XX, mod_cpp(n - K + i, N), K, N) * XX(mod_cpp(n + K + i, N)) / K;
    //std::cout << "iteration " << i << ": " << r_sum << std::endl;
  }
  
  double XX_Kn_val = (-Wn_Odd_Cpp(XX, mod_cpp(n - (2 * K), N), K, N)) * Wn_Odd_Cpp(XX, mod_cpp(n - K, N), K, N) + r_sum;
  
  return XX_Kn_val;
  
}



double RHS(const arma::colvec& XX, const int& j, const int& K, const int& F){
  
  if(K % 2 == 0){
    return XX_Kn_Even_Cpp(XX, j, K) - XX(j) + F;
  } else {
    return XX_Kn_Odd_Cpp(XX, j, K) - XX(j) + F;
  }
}



//' A single iteration of the Lorenz 04 model
//'
//' @param XX initial state
//' @param F_Lor F from the Lorenz 04 model
//' @param K_Lor K from the Lorenz 04 model
//' @param dt time step
//' @param M the number of steps that make a time step
//' @param order the order of RK method to use
//' @return dx the increment to reach the new state
//' @export
// [[Rcpp::export]]
arma::colvec DeltaLorenz04M2Cpp(const arma::colvec X0, const int & F_Lor, const int& K,
                                const double& dt, const int& M, const int & order){

  // Initialize vector of zeros that will "update" XX
  arma::colvec XX = X0;
  int N_Lor = XX.n_rows;
  arma::colvec dx(N_Lor);
  dx.fill(0);
  double F = F_Lor;
  
  // Nested for loops to update XX
  for(int i = 0; i < M ; ++i){
    for(int j = 0; j < N_Lor; ++j){
      
      if( order==4 ){
        double k1 = RHS(XX, j, K, F);
        double k2 = RHS(XX + 0.5*dt*k1, j, K, F);
        double k3 = RHS(XX + 0.5*dt*k2, j, K, F);
        double k4 = RHS(XX + dt*k3, j, K, F);
          
        dx(j) =  (k1 + 2*k2 + 2*k3 + k4) * dt/6;  
      } else if( order==1 ){
        dx(j) = RHS(XX, j, K, F)*dt;
      }
        
    } 
    XX = XX + dx;
  }
  
  return XX - X0;
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
//' @return Xiter values of all iterations
//' @export
// [[Rcpp::export]]
arma::mat Lorenz04M2SimCpp(const arma::colvec& Xinit, const int& F_Lor, const int& K_Lor,
                        const double& dt, const int& M, const int& iter, const int& burn, const int & order)
{
  // Get N_Lor from Xinit
  int N_Lor = Xinit.n_rows;
  // Create Matrix to store generated data
  arma::mat Xiter(N_Lor, iter);
  
  // Create vector to allow for burn-in iterations
  arma::colvec newXburn = Xinit;
  
  int burnin = burn;

  if( burnin > 0 ){
    for(int i = 0; i < burnin; ++i){
      newXburn = newXburn + DeltaLorenz04M2Cpp(newXburn, F_Lor, K_Lor, dt, M, order);
    }
  } 
  
  Xiter.col(0) = newXburn;
  
  for(int i = 1; i < iter; ++i){
    Xiter.col(i) = Xiter.col(i-1) + DeltaLorenz04M2Cpp(Xiter.col(i - 1), F_Lor, K_Lor, dt, M, order);
  }
  
  
  return Xiter;
}



