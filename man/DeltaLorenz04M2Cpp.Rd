% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{DeltaLorenz04M2Cpp}
\alias{DeltaLorenz04M2Cpp}
\title{A single iteration of the Lorenz 04 model}
\usage{
DeltaLorenz04M2Cpp(Z0, F, K, dt, M, vectorAlgo)
}
\arguments{
\item{Z0}{initial state}

\item{F}{F from the Lorenz 04 model}

\item{K}{K from the Lorenz 04 model}

\item{dt}{time step}

\item{M}{the number of steps that make a time step}

\item{vectorAlgo}{whether to run the vector version of the algorithm}
}
\value{
dx the increment to reach the new state
}
\description{
A single iteration of the Lorenz 04 model
}
