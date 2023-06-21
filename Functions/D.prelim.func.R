####################################
#  PART 0: PRELIMINARY FUNCTIONS   #
####################################

# 0.1. A function generating precision matrix of delta
cov.mat.model.D <- function(p, M, deltas=seq(1,30,0.5)){
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
    if(min(eig.vals) > 0)
      return(Theta)
  }
  
  stop("Model D Theta: deltas not large enough")
}