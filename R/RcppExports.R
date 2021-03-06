# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' A single iteration of the Lorenz 04 model
#'
#' @param Z0 initial state
#' @param F F from the Lorenz 04 model
#' @param K K from the Lorenz 04 model
#' @param dt time step
#' @param M the number of steps that make a time step
#' @param vectorAlgo whether to run the vector version of the algorithm
#' @return dx the increment to reach the new state
#' @export
DeltaLorenz04M2Cpp <- function(Z0, F, K, dt, M, vectorAlgo) {
    .Call('_Lorenz04_DeltaLorenz04M2Cpp', PACKAGE = 'Lorenz04', Z0, F, K, dt, M, vectorAlgo)
}

#' Simulate from the Lorenz 04 model
#'
#' @param Xinit initial state
#' @param F_Lor F from the Lorenz 04 model
#' @param K_Lor K from the Lorenz 04 model
#' @param dt time step
#' @param M the number of steps that make a time step
#' @param iter number of iterations
#' @param burn how many steps to discard
#' @param vectorAlgo whether to run the vector version of the algorithm
#' @return Xiter values of all iterations
#' @export
Lorenz04M2SimCpp <- function(Xinit, F_Lor, K_Lor, dt, M, iter, burn, vectorAlgo) {
    .Call('_Lorenz04_Lorenz04M2SimCpp', PACKAGE = 'Lorenz04', Xinit, F_Lor, K_Lor, dt, M, iter, burn, vectorAlgo)
}

#' Gradient for the Lorenz model
#'
#' @param X initial state
#' @param K K from the Lorenz 04 model
#' @param M the number of steps that make a time step
#' @param dt time step
#' @param F F from the Lorenz 04 model
#' @export
exactGradient <- function(X, K, M, dt, F) {
    .Call('_Lorenz04_exactGradient', PACKAGE = 'Lorenz04', X, K, M, dt, F)
}

#' Calculate the exact gradient in the Lorenz model
#'
#' @param X the vector at which to evaluate the gradient
#' @param K parameter from the Lorenz model
#' @return grad the matrix with the gradient
#' @export
exactGrad <- function(X, K) {
    .Call('_Lorenz04_exactGrad', PACKAGE = 'Lorenz04', X, K)
}

#' right hand side of the Lorenz model in vector form
#'
#' @param X current state
#' @param K parameter from the Lorenz 04 model
#' @param F parameter from the Lorenz 04 model
#' @export
vectorRHS <- function(X, K, F) {
    .Call('_Lorenz04_vectorRHS', PACKAGE = 'Lorenz04', X, K, F)
}

