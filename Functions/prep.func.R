# Packages
library(mvtnorm)

# Vectorize Y: from n*M to nm*1
Y.vectorize <- function(Y){
  n <- nrow(Y)
  M <- ncol(Y)
  res <- numeric(0)
  for(i in 1:n){
    res <- rbind(res, matrix(Y[i,],M,1))
  }
  res <- as.vector(res)
  return(res)
}

# Kronecker product X with I to get Z: nM*pM^2.
# Then combine a column of 1 as intercept
X.kprod <- function(X,M){
  Z <- kronecker(X, diag(M))
  return(Z)
}

# A.array to beta
col.vectorize <- function(A.array){
  R <- nrow(A.array)
  C <- ncol(A.array)
  res <- numeric(0)
  for(i in 1:C){
    res <- cbind(res, matrix(A.array[,i],R,1))
  }
  res <- as.vector(res)
  return(res)
}