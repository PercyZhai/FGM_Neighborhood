####################################
#  PART 0: PRELIMINARY FUNCTIONS   #
####################################

# 0.1. A function generating precision matrix of delta
cov.mat.model3 <- function(p, M, deltas=seq(1,30,0.5)){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN, SYMMETRY!
  Theta.off.diag <- matrix(0, nrow=p*M, ncol=p*M)
  for(i in 1:p){
    for(j in 1:p){
      if(i<j & runif(1)<0.1)
        Theta.off.diag[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0.5*diag(M)
      else
        Theta.off.diag[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
    }
  }
  # symmetrize
  for(i in 1:p){
    for(j in 1:p){
      if(i>j)
        Theta.off.diag[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- Theta.off.diag[((j-1)*M+1):(j*M), ((i-1)*M+1):(i*M)]
    }
  }
  for(delta in deltas){
    Theta.on.diag <- delta * diag(p*M)
    Theta <- Theta.on.diag + Theta.off.diag
    eig.vals <- eigen(Theta)$values
    if(min(eig.vals) > 0) return(Theta)
  }
  
  stop("Model 3 Theta: deltas not large enough")
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