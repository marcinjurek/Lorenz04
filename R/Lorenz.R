# Rfile to export Lorenz functions from cpp without allowing access to those functions
# And with wrapping to check compatability

#' Title
#'
#' @param X0 The vector of points you want to pass through a Lorenz Chaos Model.
#' @param F_Lor A Forcing constant for the Chaos Model. Recommended values are
#' @param K_Lor 
#' @param dt 
#' @param M 
#'
#' @return The function returns a new vector that represents a single step in a Lorenz chaos
#' model as described in the paper above
#' @export
#'
#' @examples
#' X0 = rnorm(960)
#' Lorenz04M2(X0, 10, 32, 0.005, 40)
#' 
Lorenz04M2 = function(X0, F_Lor, K_Lor, dt = 0.005, M = 40){
  N_Lor = length(X0)
  
  if(!is.vector(X0)){ stop("X0 should be a vector")}
  
  Xout = Lorenz04M2Cpp(X0, F_Lor, K_Lor, dt, M, N_Lor)
  return(Xout)
}

Lorenz04M2Sim = function(X0, F_Lor, K_Lor, dt, M, iter = 500, burn = 100){
  N_Lor = length(N_Lor)
  
  if(!is.vector(X0)){ stop("X0 should be a vector")}
  
  Xout = Lorenz04M2SimCpp(X0, F_Lor, K_Lor, dt, M, N_Lor, iter, burn)
  return(Xout)
}

