#######################################
#  Turning time series into FPC scores
#######################################

library(fda)

load("GitHub/GroupLasso/Real data analysis/time_series.RData")

fpc.real <- function(h, M=5){
  # Number of basis functions
  n <- dim(h)[1]
  p <- dim(h)[2]
  tau <- dim(h)[3]
  
  obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta
  
  fpc.score <- numeric(0)
  for(j in 1:p){
    obs.val.matrix <- matrix(0, nrow=tau, ncol=n)
    for (i in c(1:n)){
      obs.val.vec <- as.vector(h[i, j, ])
      obs.val.matrix[, i] <- obs.val.vec
    }
    bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
    fd.object.array <- Data2fd(argvals=obs.time, y=obs.val.matrix, basisobj=bspline.basis)
    fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
  }
  
  return(fpc.score)
}

p <- 176; M <- 5
fpc.score.1 <- fpc.real(h.1)
fpc.score.2 <- fpc.real(h.2)

save(p,M,fpc.score.1, fpc.score.2, file="Github/GroupLasso/Real data analysis/fpc_score.RData")

