## Functions for update steps
library(matrixStats)
library(invgamma)
library(truncnorm)
library(MASS)
library(lattice)
library(inline)
library(parallel)
library(igraph)
library(reshape2)
library(Matrix)
library(adaptMCMC)

# Base update function
g.mult=function(x.samp, y.samp, K){
  x.samp + K %*% (y.samp - H.mat %*% x.samp)
}

# Simple EnKF update function
simple.update=function(x.mult.prior, y.mult.i, sigma.mat, H.mat, tau.mat){
  x.mult.prior = as.matrix(x.mult.prior)
  y.mult.i = as.matrix(y.mult.i)
  temp = H.mat %*% sigma.mat %*% t(H.mat) + tau.mat
  K = sigma.mat %*% t(H.mat) %*% solve(temp)
  x.mult.post = g.mult(x.mult.prior, y.mult.i, K)
  return(x.mult.post)
}

# EnKF update with cholesky decomps
chol.update=function(x.mult.prior,y.mult.i,chol.mat, H.mat, tau.mat){
  new.sigma.chol = crossprod(chol.mat)
  
  temp.chol = new.sigma.chol + crossprod(H.mat, solve(tau.mat, H.mat))
  post.chol.mat = chol(temp.chol)
  
  temp.step = new.sigma.chol %*% x.mult.prior + crossprod(H.mat, solve(tau.mat, y.mult.i))
  step.chol = solve(t(post.chol.mat), temp.step)
  x.mult.post.chol = solve(post.chol.mat, step.chol)
  return(x.mult.post.chol)
}

# # EnKF update with Vecchia
# vec.update=function(x.mult.prior,y.mult.i,m=1,S){
#   # Vecchia Prior EnKF update
#   ## set up priors from Brian's code
#   n=nrow(x.mult.prior)
#   
#   NNarray=find_ordered_nn(S,m)
#   NNarray=as.matrix(NNarray[,-1])
#   
#   initt <- init_txx2(t(x.mult.prior), NNarray)
#   
#   thets=optim(c(1,0,-1), loglikeli, initt=initt, NNarray=NNarray, N = N, method="L-BFGS-B",lower=-5, upper=4)$par
#   thetps=thetas_to_priors(thets,n^2,m)
#   alpha=thetps[[1]]
#   beta=thetps[[2]]
#   gamma=thetps[[3]]
#   
#   ## Generate posterior values from sourced Brian's code
#   # This sourced function generates the posterior parameters
#   posts=get_posts(t(x.mult.prior),alpha,beta,gamma,NNarray)
# 
#   # This function use posterior parameters to sample from posteriors
#   # Returns an upper triangular matrix
#   vec.chol.mat=as.matrix(samp_posts(posts, NNarray))
# 
#   
#   temp.vec=vec.chol.mat%*%t(vec.chol.mat)+t(H.mat)%*%tau.mat.inv%*%H.mat
#   post.vec.mat=chol(temp.vec)
#   new.sigma.vec=vec.chol.mat%*%t(vec.chol.mat)
#   temp=new.sigma.vec%*%x.mult.prior+t(H.mat)%*%tau.mat.inv%*%y.mult.i
#   step.vec=solve(t(post.vec.mat),temp)
#   x.mult.post.vec=solve(post.vec.mat,step.vec)
# }
# 
