# Rfile to export Lorenz functions from cpp without allowing access to those functions
# And with wrapping to check compatability

Lorenz04M2 = function(X0, F_Lor, K_Lor, dt, M){
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

