library(fda)

# 1.1. Fourier basis function
fourier.basis <- function(t, M, l, omega=1){
  # Input:
  #   t, a vector of time points
  #   M bases in total
  #   l th basis function value is being calculated
  #   omega, the frequency. Default is 1.
  # Output:
  #   a vector as long as t, of fourier basis function values
  n <- length(t)
  b <- numeric(n)
  for (i in 1:n){ # basis function defined on [0,1]
    if (t[i]>=0 & t[i]<=1){
      if(l == 1) b[i] <- 1
      else if(l %% 2 == 0) b[i] <- sqrt(2) * cos(2 * pi * omega * (l/2) * t[i])
      else b[i] <- sqrt(2) * sin(2 * pi * omega * ((l-1)/2) * t[i])
    } else {
      b[i] <- 0
    }
  }
  return(b)
}

# 1.2. Observed basis function values given a vector of time points
obs.fourier.bases <- function(obs.time, M, omega=1){
  # Input:
  #   obs.time, the observation grid of length tau
  #   M, number of basis functions
  #   omega, the frequency of the fourier functions. Default: 1.
  # Output:
  #   tau * M matrix of observed function values
  tau <- length(obs.time)
  ObservMat <- matrix(0, nrow=n, ncol=M)
  for (k in 1:tau){
    for (l in 1:M){
      ObservMat[k, l] <- fourier.basis(obs.time[k], M, l, omega)
    }
  }
  return(ObservMat)
}

# 2.1.
fda.fourier.mat <- function(obs.time, M, period=1){
  # period needs to be <=1 to make the problem meaningful
  basis.fourier <- create.fourier.basis(rangeval = c(0,period), nbasis=M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj=basis.fourier)
  return(ObservMat)
}

fda.bspline.mat <- function(obs.time, M){
  # obs.time has to be from 0 to 1
  basis.bspline <- create.bspline.basis(rangeval = c(0,1), nbasis=M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj=basis.bspline)
  return(ObservMat)
}

fda.exp.mat <- function(obs.time, M){
  basis.exp <- create.exponential.basis(nbasis=M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj=basis.exp)
  return(ObservMat)
}

fda.monomial.mat <- function(obs.time, M){
  basis.monomial <- create.monomial.basis(nbasis = M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj = basis.monomial)
  return(ObservMat)
}

fda.power.mat <- function(obs.time, M){
  basis.power <- create.power.basis(nbasis = M)
  ObservMat <- eval.basis(evalarg=obs.time, basisobj = basis.power)
  return(ObservMat)
}
