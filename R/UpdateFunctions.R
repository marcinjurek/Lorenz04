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

# EnKF update with tapered cov matrix
tap.update = function(x.mult.prior, y.mult.i, sigma.mat, H.mat, tau.mat, taper.radius){
  taper = Wendland(dist.mat, taper.radius, 1, 1)
  taper.sigma.mat = sigma.mat * taper
  
  simple.update(x.mult.prior, y.mult.i, taper.sigma.mat, H.mat, tau.mat)
}