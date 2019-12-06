# Vecchia subfunctions - NOT EXPORTED


# Given Data, thetas, Nearest neighbor array, and it returns posterior 
get_posts <- function(x, a, b, g, NNarray) {
  # Get dims
  n2 <- ncol(x)
  N <- nrow(x)
  m <- min(ncol(g), ncol(NNarray))
  
  # Initialize
  a_post <- rep(0,n2)
  b_post <- rep(0,n2)
  muhat_post <- matrix(NA,nr=n2,nc=m)
  G_post <- array(NA, dim=c(m,m,n2))
  
  # Calc a_post
  a_post <- a + N/2
  
  # First b_post
  b_post[1] <- b[1] + t(x[,1]%*%x[,1])/2
  
  # Loop to calc Gamma post and b_post
  for(i in 2:n2){
    gind <- na.omit(NNarray[i,1:m])
    nn <- length(gind)
    xi <- -x[,gind]
    yi <- x[,i]
    Ginv <- t(xi)%*%xi + diag(g[i,1:nn]^(-1), nrow = nn)
    Ginv_chol <- chol(Ginv)
    #G <- ginv(Ginv)
    #muhat <- G%*%t(xi)%*%yi
    muhat <- tryCatch({solve(Ginv_chol, solve(t(Ginv_chol),crossprod(xi,yi)))},
                      error=function(e) {
                        ginv(Ginv)%*%crossprod(xi,yi)
                      })
    #G_post[1:nn,1:nn,i] <- G
    muhat_post[i,1:nn] <- muhat
    G_post[1:nn,1:nn,i] <- Ginv_chol
    muhat_post[i,1:nn] <- muhat
    b_post[i] <- b[i] + (t(yi)%*%yi - t(muhat)%*%Ginv %*% muhat)/2
  }
  return(list(a_post,b_post,muhat_post,G_post))
}


# Given posteriors, nearest neighbor info and returns samples from posteriors
samp_posts <- function(posts, NNarray, m = m) {
  # Get dims
  n2 <- nrow(NNarray)
  
  # Calculate d
  d <- (1/sqrt(posts[[2]][1]))*exp(lgamma((2*posts[[1]][1]+1)/2) - lgamma(posts[[1]][1]))
  
  # Initialize uhat
  uhat <- sparseMatrix(i=1,j=1,x=d, dims=c(n2,n2),triangular=TRUE)
  
  # Loop to calculate uhat
  for(i in 2:n2) {
    gind <- na.omit(NNarray[i, 1:m])
    
    nn <- length(gind)

    uhat[i,i] <- (1/sqrt(posts[[2]][i]))*exp(lgamma((2*posts[[1]][i]+1)/2) - lgamma(posts[[1]][i]))

    uhat[gind,i] <- posts[[3]][i,1:nn]*uhat[i, i]
  }
  return(uhat)
}



init_txx2 <- function(data, NNarray) {
  txx_txy_tyy <- list()
  txx_txy_tyy[[1]] <- list(NA,NA,crossprod(data[,1]))
  for(z in 2:nrow(NNarray)){
    gind <- na.omit(NNarray[z,])
    txx <- crossprod(-data[,gind])
    txy <- crossprod(-data[,gind],data[,z])
    txx_txy_tyy[[z]] <- list(txx,txy,crossprod(data[,z]))
  }
  return(txx_txy_tyy)
}

loglikeli <- function(x, datum, NNarray, eps, m = NULL){
  n2 <- nrow(NNarray)
  ap <- 6 + N / 2
  pr <- thetas_to_priors(x, n2, eps = eps, m)
  if(is.null(m)){
    m <- min(ncol(pr[[3]]), ncol(NNarray))
  }
  if(m < 2){m <- 2}
  b <- pr[[2]]; g <- pr[[3]];
  sums <- 6*log(b[1]) - ap*log(b[1] + crossprod(datum[,1])/2)
  for(i in 2:n2) {
    gind <- na.omit(NNarray[i,1:m])
    nn <- length(gind)
    # browser()
    xi <- -datum[,gind]
    yi <- datum[,i]
    Ginv <- tryCatch({crossprod(xi) + diag(g[i,1:nn]^(-1),nrow=nn)},
                     error=function(e) {show(dim(crossprod(xi))); show(g[i,1:nn]);show(m);show(nn);show(length(na.omit(NNarray[i,1:m])));show(i)})
    # Ginv <- crossprod(xi) + diag(g[i,1:nn]^(-1),nrow=nn)
    Ginv_chol <- chol(Ginv)
    #G <- ginv(Ginv)
    #muhat <- G%*%t(xi)%*%yi
    #muhat <- solve(Ginv_chol, solve(t(Ginv_chol),initt[[i]][[2]]))
    muhat <- tryCatch({solve(Ginv_chol, solve(t(Ginv_chol),crossprod(xi,yi)))},
                      error=function(e) {
                        ginv(Ginv)%*%crossprod(xi,yi)
                      })
    #G_post[1:nn,1:nn,i] <- G
    b_post <- b[i] + (crossprod(yi) - t(muhat)%*%Ginv %*% muhat)/2
    ldet <- -0.5*(2*sum(log(na.omit(diag(Ginv_chol)))) + (sum(log(g[i,1:nn]))))
    lb <- 6*log(b[i]) - ap*log(b_post)
    sums <- sums + ldet+lb
  }
  return(c(-sums))
}


thetas_to_priors <- function(thetas, n2, eps, m = NULL) {
  b <- 5 * exp(thetas[[1]])*(1 - exp(-exp(thetas[[2]])/sqrt(0:(n2 - 1))))
  a <- rep(6, n2)
  theta_3 <- exp(-exp(thetas[[3]]) * (1:30))
  if(is.null(m)){
    m <- which(theta_3 < eps)[1] - 1
  }
  if(is.na(m) | m < 2){m <- 2}
  g <- matrix(theta_3[1:m], nc = m, nr = n2, byrow = T)
  g <- g / (b / (a - 1))
  return(list(a, b, g))
}