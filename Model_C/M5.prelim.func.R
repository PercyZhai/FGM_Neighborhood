####################################
#  PART 0: PRELIMINARY FUNCTIONS   #
####################################

# 0.1. A function generating precision matrix of delta
cov.mat.model5 <- function(p,M){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN
  Theta <- matrix(nrow=p*M, ncol=p*M)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M)
      else if(abs(i-j)==1) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.4
      else if(abs(i-j)==2) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.2
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
  ObservMat <- matrix(0, nrow=tau, ncol=M)
  for (k in 1:tau){
    for (l in 1:M){
      ObservMat[k, l] <- fourier.basis(obs.time[k], M, l)
    }
  }
  return(ObservMat)
}

library(poweRlaw)

# Step 1. Generate edge set E
# output: an adjacency matrix G
# using power law distribution
edge.gen <- function(pi=0.05, p, alpha=2, seed){
  
  set.seed(seed)
  G <- diag(p)
  num.edg <- 0
  
  # G is symmetric, so select the neighbors in the upper and lower triangle
  for (i in 1:(p-1)){
    # number of neighbors: follows power law distribution
    d <- rpldis(1, xmin=1, alpha=alpha)
    if (d > p-i) d <- p-i # maximum neighbors selectable: p-i
    # select the exact neighbors of this node
    neigh <- sample((i+1):p, size=d, replace=F)
    
    for (k in 1:d){
      G[i, neigh[k]] <- 1
      G[neigh[k], i] <- 1
    }
    
    # count the number of edges up till now
    num.edg <- num.edg + d
    if (num.edg > p*(p-1)/2 * pi){
      return(G)
    }
  }
  
  return(G)
}


# Step 2. Generate E1,...,EM from E
# input: adjacency matrix G
# output: adjacency matrices G_l, l=1,...,M
edge.sets <- function(G, tau=0, M=20, seed){
  p <- nrow(G)
  G.list <- list()
  
  set.seed(seed)
  
  # 1
  # Generating Ec's adjacency matrix Gc
  Gc <- diag(p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j & G[i,j]==1 & runif(1)<tau) Gc[i,j] <- 1
    }
  }
  # 2
  for(l in 1:M) G.list[[l]] <- Gc
  # 3
  l <- 1; B <- 1
  # 4-10
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(G[i,j]==1 & Gc[i,j]!=1){
        G.list[[l]][i,j] <- 1
        G.list[[l]][j,i] <- 1
        l <- l+1
        if(l>B){
          l <- 1
          B <- (B+1) %% M
        }
      }
    }
  }
  return(G.list)
}



tilde.to.omega <- function(A){
  p <- nrow(A)
  # step 1. rescale by row: divided by 2-norm
  for(j in 1:p){
    A[,j] <- A[,j] / sqrt(sum((A[,j])^2))
  }
  
  # step 2. symmetrize: average with transpose
  A <- (A + t(A))/2
  
  # step 3. diagonal entries to one
  for(i in 1:p){
    A[i,i] <- 1
  }
  return(A)
}


# Function to generate Omega1 to OmegaM from E1,...,EM, given G.list[[1]] to[[M]]
Omega.gen <- function(G.list, seed){
  set.seed(seed)
  M <- length(G.list)
  p <- nrow(G.list[[1]])
  Omega.tilde.list <- list()
  for(l in 1:M){
    Omega.tilde.list[[l]] <- diag(p)
    for(j in 1:(p-1)){
      for(i in (j+1):p){
        if(G.list[[l]][i,j]==1){
          u.D <- runif(1, 1/3, 2/3)
          if(runif(1)<0.5) u.D <- -u.D
          Omega.tilde.list[[l]][i,j] <- u.D
        }
      }
    }
  }
  # Now we have Omega.tilde.list.
  Omega.list <- list()
  for(l in 1:M){
    Omega.list[[l]] <- tilde.to.omega(Omega.tilde.list[[l]])
  }
  return(Omega.list)
}




# Function to return Sigma of models 4 and 5 given a list of omega
Sigma.mod.4 <- function(Omega.list){
  M <- length(Omega.list)
  p <- nrow(Omega.list[[1]])
  Sigma <- matrix(0, p*M, p*M)
  
  for(l in 1:M){
    Sigma[(l-1)*p + 1:p, (l-1)*p + 1:p] <- 3 * l^(-1.8) * solve(Omega.list[[l]])
  }
  return(Sigma)
}

Sigma.mod.5 <- function(Omega.list){
  Sigma.ps <- Sigma.mod.4(Omega.list)
  M <- length(Omega.list)
  p <- nrow(Omega.list[[1]])
  Omega <- matrix(0, p*M, p*M)
  for(l in 1:M){
    if(l==M) Omega[(l-1)*p + 1:p, (l-1)*p + 1:p] <- Omega.list[[l]]
    else{
      Omega[(l-1)*p + 1:p, (l-1)*p + 1:p] <- Omega.list[[l]]
      Omega.l.star <- Omega.list[[l]] - diag(diag(Omega.list[[l]]))
      Omega.lp1.star <- Omega.list[[l+1]] - diag(diag(Omega.list[[l+1]]))
      Omega.off.diag <- (Omega.l.star + Omega.lp1.star)/2
      Omega[(l)*p + 1:p, (l-1)*p + 1:p] <- Omega.off.diag
      Omega[(l-1)*p + 1:p, (l)*p + 1:p] <- Omega.off.diag
    }
  }
  S.diag.rt <- diag(sqrt(diag(Sigma.ps)))
  O.diag.invrt <- solve(diag(sqrt(diag(Omega))))
  Sigma.nps <- S.diag.rt %*% solve(O.diag.invrt %*% Omega %*% O.diag.invrt) %*% S.diag.rt
  return(Sigma.nps)
}