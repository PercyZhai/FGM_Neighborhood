# Qiao's method
library(matrixcalc)

blockwise_Frob <- function(ma, M){
  # ma: matrix
  # block size
  p <- dim(ma)[1]
  if (p %% M != 0){
    stop("dimension of matrix cannot be divided by block size")
  }
  n.b <- p / M
  result <- matrix(0, nrow=n.b, ncol=n.b)
  for (i in 1:n.b){
    for (j in 1:n.b){
      result[i, j] <- norm(ma[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)], "F")
    }
  }
  return(result)
}

object_func <- function(S, Theta, gamma, M){
  
  part1 <- -log(det(Theta))
  part2 <- sum(diag(S %*% Theta))
  
  block_frob <- blockwise_Frob(Theta, M)
  part3 <- gamma * (sum(block_frob) - sum(diag(block_frob)))
  
  return(part1+part2+part3)
}

ProxAlg_FGM <- function(S, p, M, gamma, Eta="Auto", n.iteration=2000){
  # Optimzation function to solve the optimization problem in FGM
  # Use Proximal Algorithm
  # 
  # Args:
  #   S: the estimated covariance matrix of Graph
  #   p: the number of vertices
  #   M: the number of principle components
  #   gamma: hyperparameter.
  #   Eta: hyperparamter. If Eta="Auto", then Eta would be set as inverse of the 
  #        product of the maximum singular value of cov.X and cov.Y
  #   n.iteration: maximum times of iteration
  #
  # Returns:
  #   A list:
  #     ThetaMathat: The estimated differential matrix between ThetaX and ThetaY
  #     blockFrob: A matrix of Frobenius norm of each block matrix
  #     Support: Support Matrix
  if(Eta == "Auto"){
    Eta <- 0.01 # this needs more careful design
  }
  
  v.obj <- c()
  
  converge.indicator <- FALSE
  
  Theta.old <- diag(1, nrow=(p*M), ncol=(p*M))
  for (t in 1:n.iteration){
    obj.old <- object_func(S, Theta.old, gamma, M)
    v.obj <- c(v.obj, obj.old)
    
    Theta <- Theta.old
    
    # calculate gradient
    grad.Mat <- S - round(solve(Theta), 10)  # round to solve asymmetric problem due to 
    # numerical issue
    A <- Theta - Eta * grad.Mat
    
    # update Theta
    for (i in 1:p){
      for (j in 1:p){
        Asub <- A[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]
        Fnorm <- norm(Asub,"F")
        if (i != j){
          if(Fnorm <= gamma*Eta){
            Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, nrow=M, ncol=M)
          } else {
            Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- round(((Fnorm - gamma*Eta) / Fnorm) * Asub, 9)
          }
        } else {
          Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- round(Asub, 9)
        }
      }
    }
    
    Theta.new <- Theta
    obj.new <- object_func(S, Theta.new, gamma, M)
    
    if(is.na(obj.new)| obj.new == -Inf){
      Theta.new <- Theta.old
      Theta.hat <- 0.5 * (Theta.new + t(Theta.new))
      Theta.frob <- blockwise_Frob(Theta.hat, M)
      
      SupMat <- matrix(0, nrow=p, ncol=p)
      SupMat[Theta.frob > 0] <- 1
      
      return (list(ThetaMathat=Theta.hat, blockFrob=Theta.frob, Support=SupMat, 
                   converge=converge.indicator, num.iter=t))
      break
    }
    
    if (abs((obj.new-obj.old)/(obj.old+1e-15)) < 1e-3){
      converge.indicator <- TRUE
      Theta.hat <- 0.5 * (Theta.new + t(Theta.new))
      Theta.frob <- blockwise_Frob(Theta.hat, M)
      
      SupMat <- matrix(0, nrow=p, ncol=p)
      SupMat[Theta.frob > 0] <- 1
      
      return (list(ThetaMathat=Theta.hat, blockFrob=Theta.frob, Support=SupMat, 
                   converge=converge.indicator, num.iter=t))
      break
    }
    Theta.old <- Theta.new
  }
  
  Theta.hat <- 0.5 * (Theta.new + t(Theta.new))
  Theta.frob <- blockwise_Frob(Theta.hat, M)
  
  SupMat <- matrix(0, nrow=p, ncol=p)
  SupMat[Theta.frob > 0] <- 1
  
  
  return (list(ThetaMathat=Theta.hat, blockFrob=Theta.frob, Support=SupMat, 
               converge=converge.indicator, num.iter=t, min.obj=obj.new))
}