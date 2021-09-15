####################################
#  PART 0: PRELIMINARY FUNCTIONS   #
####################################

# Generating M*M Tridiagonal
tridiag <- function(M){
  result <- diag(M)
  for(i in 1:M){
    for(j in 1:M){
      if(abs(i-j)==1) result[i,j] <- 0.5
    }
  }
  return(result)
}

# 0.1. A function generating precision matrix of delta
cov.mat.model6 <- function(p,M){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN
  Theta <- matrix(nrow=p*M, ncol=p*M)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M)
      else if(abs(i-j)==1) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.4
      else if(abs(i-j)==2) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.2
      else Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
    }
  }
  return(Theta)
}

# 0.2. Fourier basis function
fourier.basis <- function(t, M, l){
  # Input:
  #   t, a vector of time points
  #   M bases in total
  #   l th basis function value is being calculated
  # Output:
  #   a vector as long as t, of fourier basis function values
  n <- length(t)
  b <- numeric(n)
  for (i in 1:n){ # basis function defined on [0,1]
    if (t[i]>=0 & t[i]<=1){
      if(l == 1) b[i] <- 1
      else if(l %% 2 == 0) b[i] <- sqrt(2) * cos(2 * pi * (l/2) * t[i])
      else b[i] <- sqrt(2) * sin(2 * pi * ((l-1)/2) * t[i])
    } else {
      b[i] <- 0
    }
  }
  return(b)
}

# 0.3. Observed basis function values given a vector of time points
obs.fourier.bases <- function(obs.time, M){
  # Input:
  #   obs.time, the observation grid of length tau
  #   M, number of basis functions
  # Output:
  #   tau * M matrix of observed function values
  tau <- length(obs.time)
  ObservMat <- matrix(0, nrow=n, ncol=M)
  for (k in 1:tau){
    for (l in 1:M){
      ObservMat[k, l] <- fourier.basis(obs.time[k], M, l)
    }
  }
  return(ObservMat)
}