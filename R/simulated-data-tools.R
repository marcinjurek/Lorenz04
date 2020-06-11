getX0 = function(N, Force, K, dt, dir = '~/HVLF/simulations-lorenz/'){
  
  fileName = paste("init_Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep="_")
  filePath = paste(dir, fileName, sep="")
  
  generateInit = !file.exists(filePath)

  if( generateInit ){
    message("Initial state for these parameters has not yet been generated. Generating now.")
    X0 = rnorm(N)  
    X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=100, burn=0)
    x0 = X1[,100]
    write(x0, file=filePath)  
  } else {
    x0 = scan(filePath, quiet = TRUE)  
  }
  if(class(x0)=='numeric'){
    x0 = matrix(x0, ncol=1)
  }
  return(x0)

}



center_operator = function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}



getLRMuCovariance = function(N, Force, dt, K){

  fileName.all = paste("~/HVLF/simulations-lorenz/Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
  X = Matrix::Matrix(scan(fileName.all, quiet = TRUE), nrow = N)
  Xbar = matrix(rowMeans(as.matrix(X)), ncol=1)
  X = center_operator(X)
  S = matrix((X %*% Matrix::t(X)) / (ncol(X) - 1), ncol=N)

  return(list(mu = Xbar, Sigma = S))

}
