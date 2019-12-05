## Vecchia Update Source File


###########################################################
##   This file contains only the Vecchia Update function ##
## and all of its necessary subfunctions. 
###########################################################


# EnKF update function with Vecchia covariance approximation
vec.update = function(x.mult.prior, y.mult.i, H.mat, tau.mat, S, eps = 1e-03, NN = 30, m = NULL){
  # Vecchia Prior EnKF update
  ## set up priors from Brian's code
  ## Get the number of rows in our x.mult.prior
  n = nrow(x.mult.prior)
  x.mult.prior = as.matrix(x.mult.prior)
  
  ## Make sure eps is within a certain range
  if(eps < 1e-06 | eps > 5e-01)
  {
    stop("eps is not an acceptable value. Might cause numerical issues")
  }
  
  ## Have to generate our Nearest Neighbors array
  
  NNarray = find_ordered_nn(S, min(NN, m))
  NNarray = as.matrix(NNarray[, -1])
  
  # Use the data and nearest neighbor info to get initial values for theta 
  # log-likelihood optimization
  # initt = init_txx2(t(x.mult.prior), NNarray)
  
  # Get maximum a posteriori values for theta
  thets = optim(c(1,0,-1), loglikeli, datum = t(x.mult.prior), NNarray = NNarray, 
                eps = eps, m = m, method="L-BFGS-B", lower=-2, upper=3)$par
  
  # Use M.a.p. theta values to generate prior values for alpha, beta, and gamma
  thetps = thetas_to_priors(thets, n^2, eps = eps, m)
  alpha = thetps[[1]]
  beta = thetps[[2]]
  gamma = thetps[[3]]
  m = min(ncol(gamma), ncol(NNarray))
  
  ## Generate posterior values from sourced Brian's code
  # This sourced function generates the posterior parameters
  posts = get_posts(t(x.mult.prior), alpha, beta, gamma, NNarray)
  
  # This function use posterior parameters to sample from posteriors
  # Returns an upper triangular matrix
  vec.chol.mat = as.matrix(samp_posts(posts, NNarray, m))
  
  # Use chol.update function to generate our update?
  
  ## Calculate update using formula described in paper
  # Calculate new sigma and t(H.mat)%*%tau.mat.inv using crossprod because they 
  # need to be used twice in the code beneath
  tau.mat.inv=solve(tau.mat)
  new.sigma.vec = tcrossprod(vec.chol.mat)
  H.tau.mat = crossprod(H.mat, tau.mat.inv)
  
  # Calculations using update formula
  temp.vec = new.sigma.vec + H.tau.mat %*% H.mat
  post.vec.mat = chol(temp.vec)
  temp = new.sigma.vec %*% x.mult.prior + H.tau.mat %*% y.mult.i
  
  # Use solve against a vector instead of just a matrix to speed up calculation
  step.vec = solve(t(post.vec.mat), temp)
  x.mult.post.vec = solve(post.vec.mat, step.vec)
  
  # Return the updated sample
  return(list("update" = x.mult.post.vec, "m" = m))
}
