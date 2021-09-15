#####################################
# Load FPC score matrix and compute adjacency matrix
#####################################

# get M, p, fpc.score.1 and 2
load('/home/percysjzhai/GroupLasso/Real data analysis/time_series.RData')
h.1 <- h.1[1:70,,]
h.2 <- h.2[1:95,,]
n.1 <- dim(h.1)[1]
n.2 <- dim(h.2)[1]
p <- dim(h.2)[2]
L <- 50
M <- 5

library(fda)

func.path <- "/home/percysjzhai/GroupLasso/Functions"
source(paste(func.path,"prep.func.R", sep="/"))
source(paste(func.path,"ADMM.R", sep="/"))
source(paste(func.path,"grplasso.iter.R", sep="/"))
source(paste(func.path,"scv.bic.R", sep="/"))
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))
source(paste(func.path,"prec.rec.R", sep="/"))
source(paste(func.path,"FPCA.score.R", sep="/"))

## code copied from adjmat.M1p50.R

A.mean.mat <- matrix(NA,p,p-1) # a matrix saving the mean of A
G.mat <- matrix(NA, p,p) # adjacency matrix


time.start <- proc.time()
# parallel computing
library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)

G.mat.1 <- foreach(j=1:p, .combine="rbind", .packages="fda") %dopar% {
  Y.1 <- FPC.score.j(h.1, j, M)$Y
  X.1 <- FPC.score.j(h.1, j, M)$X
  
  lambda.max <- lambda.sup(X.1,Y.1)
  lambdas <- exp(seq(log(lambda.max), log(1), length.out=L)) # This is key: Don't let lambda be too small
  
  # ADMM group lasso
  result <- group.lasso.SCV.ADMM.oneloop(X.1, Y.1, lambdas)
  
  A.list <- result$A.opt
  best.lambda <- result$lambda
  
  # plot SCV.AIC vs lambda
  aic.plot <- plot(lambdas, result$error.path, type="l", main="SCV.AIC vs lambda")
  
  # Set a threshold to recognize neighbors
  A.mean <- rep(0, p-1)
  A.nonzero <- rep(0, p-1)
  for(juliet in 1:(p-1)){
    A.mean[juliet] <- mean(A.list[[juliet]])
    A.nonzero[juliet] <- (A.mean[juliet]!=0)
  }
  if(sum(A.nonzero)==0){ # all A are 0
    threshold <- 0
  }else{
    threshold <- mean(abs(A.mean[which(A.nonzero!=0)])) * 0.05
  }
  
  # Neighbor recognization
  V.j <- rep(0, p)
  V.j[j] <- 1
  for(juliet in 1:(p-1)){
    if(juliet<j){
      if(abs(A.mean[juliet]) > threshold) V.j[juliet] <- 1
    }else{
      if(abs(A.mean[juliet]) > threshold) V.j[juliet + 1] <- 1
    }
  }
  
  V.j
}

save(G.mat.1, file="/home/percysjzhai/GroupLasso/Real data analysis/G1.RData")

G.mat.2 <- foreach(j=1:p, .combine="rbind", .packages="fda") %dopar% {
  Y.2 <- FPC.score.j(h.2, j, M)$Y
  X.2 <- FPC.score.j(h.2, j, M)$X
  
  lambda.max <- lambda.sup(X.2,Y.2)
  lambdas <- exp(seq(log(lambda.max), log(1), length.out=L)) # This is key: Don't let lambda be too small
  
  # ADMM group lasso
  result <- group.lasso.SCV.ADMM.oneloop(X.2, Y.2, lambdas)
  
  A.list <- result$A.opt
  best.lambda <- result$lambda
  
  # plot SCV.AIC vs lambda
  aic.plot <- plot(lambdas, result$error.path, type="l", main="SCV.AIC vs lambda")
  
  # Set a threshold to recognize neighbors
  A.mean <- rep(0, p-1)
  A.nonzero <- rep(0, p-1)
  for(juliet in 1:(p-1)){
    A.mean[juliet] <- mean(A.list[[juliet]])
    A.nonzero[juliet] <- (A.mean[juliet]!=0)
  }
  if(sum(A.nonzero)==0){ # all A are 0
    threshold <- 0
  }else{
    threshold <- mean(abs(A.mean[which(A.nonzero!=0)])) * 0.05
  }
  
  # Neighbor recognization
  V.j <- rep(0, p)
  V.j[j] <- 1
  for(juliet in 1:(p-1)){
    if(juliet<j){
      if(abs(A.mean[juliet]) > threshold) V.j[juliet] <- 1
    }else{
      if(abs(A.mean[juliet]) > threshold) V.j[juliet + 1] <- 1
    }
  }
  
  V.j
}

save(G.mat.2, file="/home/percysjzhai/GroupLasso/Real data analysis/G2.RData")

stopCluster(cl)


time.end <- proc.time()
time.run <- (time.end - time.start)[3]
print(time.run)