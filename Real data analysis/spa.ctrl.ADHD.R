#####################################
# Load FPC score matrix and compute adjacency matrix
# Sparsity control by alteration of threshold
#####################################

func.path <- ""
data.path <- ""     # where the saved observed data comes from
save.path <- ""     # where the final adjacency matrices are saved
runtime.path <- ""  # where runtime is saved

library(fda)
library(matrixcalc)
library(MASS)
library(mvtnorm)

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"A.prelim.func.R", sep="/")) # For Model A generation
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))   # FGLasso

# parameters
M <- 5 # number of basis used for neighborhood selection
tol.abs <- 1e-4  # Tolerance (absolute) in ADMM
tol.rel <- 1e-4  # Tolerance (relative) in ADMM
L <- 100 # number of lambdas
K <- 5 # number of folds of SCV
t.vec <- c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3, 5, 8, 10, 15, 20) # vector of t, where threshold epsilon = t * lambda
len.t <- length(t.vec)
spar.ctrl <- 0.02

time.start <- proc.time()


####################################
#   PART 1: READ IN OBSERVATION    #
####################################

#    Read in Observation h_ijk for both ADHD group and control group

# get M, p, fpc.score.1 and 2
load(paste(data.path, "time.series.ADHD.RData", sep="/"))
n.1 <- dim(h.1)[1] #ADHD group
n.2 <- dim(h.2)[1] #control group
p <- dim(h.2)[2]
tau <- dim(h.2)[3]

####################################
#     PART 2: GAIN FPC SCORE       #
####################################

# For the use of gX group
time.start.fpc <- proc.time()

fpc.score.1 <- numeric(0)
fpc.score.2 <- numeric(0)
obs.time <- seq(1/tau, 1, 1/tau)
for(j in 1:p){
  obs.val.matrix.1 <- matrix(0, nrow=tau, ncol=n.1)
  obs.val.matrix.2 <- matrix(0, nrow=tau, ncol=n.2)
  for (i in c(1:n.1)){
    obs.val.vec.1 <- as.vector(h.1[i, j, ])
    obs.val.matrix.1[, i] <- obs.val.vec.1
  }
  for (i in c(1:n.2)){
    obs.val.vec.2 <- as.vector(h.2[i, j, ])
    obs.val.matrix.2[, i] <- obs.val.vec.2
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
  # Construct a functional data object from the observation
  # bspline basis is the default setting of this function
  # It does not mean that the basis function is bspline!
  fd.object.array.1 <- Data2fd(argvals=obs.time, y=obs.val.matrix.1, basisobj=bspline.basis)
  fd.object.array.2 <- Data2fd(argvals=obs.time, y=obs.val.matrix.2, basisobj=bspline.basis)
  # FPCA process
  fpc.score.1 <- cbind(fpc.score.1, pca.fd(fd.object.array.1, nharm=M)$scores)
  fpc.score.2 <- cbind(fpc.score.2, pca.fd(fd.object.array.2, nharm=M)$scores)
}

time.end.fpc <- proc.time()
runtime.fpc <- (time.end.fpc - time.start.fpc)[3]


#######################################
##                                   ##
##      START PARALLEL COMPUTING     ##
##                                   ##
#######################################
library(doParallel)

# cores <- detectCores()

# Use the following line to show ADMM iterations for debugging
cl <- makeCluster(28, outfile="")

# Use the following line for faster direct computation
#cl <- makeCluster(64)

registerDoParallel(cl)


#######################################
#   PART 3: CACHE FOR MULTIPLE RUNS   #
#######################################

# 3.1. Get Ruiz scaling
# A.Y.1 <- matrix(NA, n.1, M)
# A.X.1 <- matrix(NA, n.1, (p-1)*M)
# 
# A.Y.2 <- matrix(NA, n.1, M)
# A.X.2 <- matrix(NA, n.1, (p-1)*M)
# 
# d.array.1 <- foreach(j = 1:p, .combine="rbind") %dopar% {
#   jth.range <- (((j-1)*M+1) : (j*M))
#   for(i in 1:n.1){
#     A.Y[i,] <- fpc.score.1[i, jth.range]
#     A.X[i,] <- fpc.score.1[i, -jth.range]
#   }
#   
#   # Ruiz 
#   d.j <- ruiz.scale(t(A.X.1) %*% A.X.1)$d
#   
#   d.j
# }
# d.array.2 <- foreach(j = 1:p, .combine="rbind") %dopar% {
#   jth.range <- (((j-1)*M+1) : (j*M))
#   for(i in 1:n.2){
#     A.Y[i,] <- fpc.score.2[i, jth.range]
#     A.X[i,] <- fpc.score.2[i, -jth.range]
#   }
#   
#   # Ruiz 
#   d.j <- ruiz.scale(t(A.X.1) %*% A.X.1)$d
#   
#   d.j
# }

d.array.1 <- matrix(1, nrow=p, ncol=(p-1)*M)
d.array.2 <- matrix(1, nrow=p, ncol=(p-1)*M)

d.1 <- list()
d.2 <- list()
norm.adj.1 <- rep(NA,p)
norm.adj.2 <- rep(NA,p)
for(k in 1:p){
  d.1[[k]] <- d.array.1[k,]
  norm.adj.1[k] <- norm(d.1[[k]],"2")
  d.2[[k]] <- d.array.2[k,]
  norm.adj.2[k] <- norm(d.2[[k]],"2")
}

################################################
##                                            ##
## PART 4. RUN ADMM-gX OVER EACH j AND lambda ##
##                                            ##
################################################

time.start.ADHD <- proc.time()

# default warm start PQU
P.def <- matrix(0, (p-1)*M, M)
Q.def <- matrix(0.1, (p-1)*M, M)
U.def <- matrix(0.01, (p-1)*M, M)

# G.mat is the optimal adjacency matrix to save
G.mat.ADHD <- matrix(NA, p, p)
G.mat.control <- matrix(NA, p, p)

A.Y <- matrix(NA, n.1, M)
A.X <- matrix(NA, n.1, (p-1)*M)

################################################
##       PART 4.1. G.mat for ADHD Group     ##
################################################

B.mat.ADHD <- foreach(j = 1:p, .combine="rbind") %dopar% {
  
  jth.range <- (((j-1)*M+1) : (j*M))
  for(i in 1:n.1){
    A.Y[i,] <- fpc.score.1[i, jth.range]
    A.X[i,] <- fpc.score.1[i, -jth.range]
  }
  
  # find maximum lambda specifically for each j
  lambdas <- seq(lambda.sup(A.X, A.Y), 0, length.out=L)
  
  # warmstart for the full dataset
  P.full <- P.def; Q.full <- Q.def; U.full <- U.def
  
  N.j <- rep(0,L)
  
  SCV.mat <- matrix(NA, L, len.t)
  B.full <- list() # a list with B.hat under each lambda.
  
  for(l in 1:L){
    lambda <- lambdas[l]
    # Step 1. Use the full dataset to estimate B.hat(lambda)
    cat(paste("Lambda ",l," of node ",j,"\n",sep=""))
    full.result <- ADMM.grplasso(A.X = A.X, A.Y = A.Y, d = d.1[[j]],
                                 lambda=lambda, rho.init=1,
                                 P.in=P.full, Q.in=Q.full, U.in=U.full,
                                 tol.abs = tol.abs, tol.rel = tol.rel,
                                 maxiter = 400)
    P.full <- full.result$P # this is B.hat(lambda)
    Q.full <- full.result$Q
    U.full <- full.result$U # will be used for warmstart for next lambda
    
    B.full[[l]] <- P.full
    
    # Process the estimated P.full into a neighborhood selection vector
    P.full.frob <- rep(0, p-1)
    for(k in 1:(p-1)){
      P.full.frob[k] <- norm(P.full[(k-1)*M + (1:M), ], "F")
    }
    
    for(ind.t in 1:len.t){ # for each epsilon
      
      # Step 2. finding N.hat.j according to each block of B.full
      threshold <- t.vec[ind.t] * lambda
      # The following corresponds to N.hat.j(lambda,epsilon)
      N.hat <- which(P.full.frob > threshold)
      # cardinality of N.hat.B for this specific lambda
      card.N.hat <- length(N.hat)
      
      slice.p <- (P.full.frob > threshold) # length p, true or false
      slice.pM <- rep(slice.p, each=M)
      
      A.X.eff <- A.X[, slice.pM] # select only the columns of A.X with selected features
      
      # Step 3. B check, related to lambda and epsilon
      B.check <- P.full
      for(k in 1:(p-1)){
        if(P.full.frob[k] <= threshold){
          # block not large enough: filter out this block
          B.check[(k-1)*M + (1:M), ] <- 0
        }
      }
      
      SCV.single <- rep(NA, K)
      
      for(k in 1:K){
        if(card.N.hat > 0){ # Circumstance A. if |N.hat| >=1
          
          # Step A4: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n.1/K) + 1, floor(k*n.1/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.X.cv <- A.X.eff[fold.k.ind, ]
          A.X.train <- A.X.eff[-fold.k.ind, ]
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step A5: calculate B.tilde from a ridge regression
          B.tilde <- solve(t(A.X.train) %*% A.X.train + 0.1*diag(card.N.hat * M)) %*% t(A.X.train) %*% A.Y.train
          
          # Step A6: evaluate SCV for this single k
          residual <- A.Y.cv - A.X.cv %*% B.tilde # residual is n.1*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat
          
        }else{ # Circumstance B. if cardinality of N.hat is 0
          # Step B4: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n.1/K) + 1, floor(k*n.1/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step B6: evaluate SCV for this single k
          residual <- A.Y.cv # residual is n.1*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          SCV.single[k] <- norm(residual,"F")^2
        }
      }# end of the k'th fold
      
      # Step 7: averaging SCV over K folds for this specific lambda and epsilon
      SCV.mat[l, ind.t] <- mean(SCV.single)
      
    } # end of all thresholds (epsilon) 1:len.t
    
  } # end of lambda 1:L
  
  # Step 8: finding the minimal SCV and its corresponding lambda/epsilon
  scv.min <- which(SCV.mat == min(SCV.mat), arr.ind=T)
  index.optimal <- scv.min[dim(scv.min)[1], ]
  l.optimal <- index.optimal[1]
  ind.t.optimal <- index.optimal[2]
  lambda.optimal <- lambdas[l.optimal]
  t.optimal <- t.vec[ind.t.optimal]
  
  SCV.optimal <- SCV.mat[index.optimal]
  B.optimal <- B.full[[l.optimal]]
  
  # Process the estimated B into a neighborhood selection vector
  B.frob <- rep(0, p-1)
  for(k in 1:(p-1)){
    B.frob[k] <- norm(B.optimal[(k-1)*M + (1:M), ], "F")
  }
  
  B.j <- rep(0,p)
  for(juliet in 1:(p-1)){
    if(juliet < j){
      B.j[juliet] <- B.frob[juliet]
    }else{
      B.j[juliet+1] <- B.frob[juliet]
    }
  }
  
  B.j
}

B.ADHD.sym <- matrix(NA, p, p)
for(j in 1:p){
  for(k in 1:p){
    B.ADHD.sym[j,k] <- min(B.mat.ADHD[j,k], B.mat.ADHD[k,j]) # equivalent to AND
  }
}

thres.vec <- seq(0, max(B.ADHD.sym)*1.01, length.out=500)

for(i in 1:500){
  threshold <- thres.vec[i]
  G.ADHD.sym <- matrix(0, p, p)
  for(j in 1:p){
    for(k in 1:p){
      G.ADHD.sym[j,k] <- as.numeric(B.ADHD.sym[j,k] > threshold)
      if(j==k) G.ADHD.sym[j,k] <- 0
    }
  }
  if(sum(G.ADHD.sym) <= spar.ctrl * p*(p-1)) break
}

time.end.ADHD <- proc.time()
runtime.ADHD <- (time.end.ADHD - time.start.ADHD)[3]

################################################
##     PART 4.2. G.mat for control Group      ##
################################################

A.Y <- matrix(NA, n.2, M)
A.X <- matrix(NA, n.2, (p-1)*M)

time.start.control <- proc.time()

B.mat.control <- foreach(j = 1:p, .combine="rbind") %dopar% {
  
  jth.range <- (((j-1)*M+1) : (j*M))
  for(i in 1:n.2){
    A.Y[i,] <- fpc.score.2[i, jth.range]
    A.X[i,] <- fpc.score.2[i, -jth.range]
  }
  
  # find maximum lambda specifically for each j
  lambdas <- seq(lambda.sup(A.X, A.Y), 0, length.out=L)
  
  # warmstart for the full dataset
  P.full <- P.def; Q.full <- Q.def; U.full <- U.def
  
  N.j <- rep(0,L)
  
  SCV.mat <- matrix(NA, L, len.t)
  B.full <- list() # a list with B.hat under each lambda.
  
  for(l in 1:L){
    lambda <- lambdas[l]
    # Step 1. Use the full dataset to estimate B.hat(lambda)
    cat(paste("Lambda ",l," of node ",j,"\n",sep=""))
    full.result <- ADMM.grplasso(A.X = A.X, A.Y = A.Y, d = d.2[[j]],
                                 lambda=lambda, rho.init=1,
                                 P.in=P.full, Q.in=Q.full, U.in=U.full,
                                 tol.abs = tol.abs, tol.rel = tol.rel,
                                 maxiter = 400)
    P.full <- full.result$P # this is B.hat(lambda)
    Q.full <- full.result$Q
    U.full <- full.result$U # will be used for warmstart for next lambda
    
    B.full[[l]] <- P.full
    
    # Process the estimated P.full into a neighborhood selection vector
    P.full.frob <- rep(0, p-1)
    for(k in 1:(p-1)){
      P.full.frob[k] <- norm(P.full[(k-1)*M + (1:M), ], "F")
    }
    
    for(ind.t in 1:len.t){ # for each epsilon
      
      # Step 2. finding N.hat.j according to each block of B.full
      threshold <- t.vec[ind.t] * lambda
      # The following corresponds to N.hat.j(lambda,epsilon)
      N.hat <- which(P.full.frob > threshold)
      # cardinality of N.hat.B for this specific lambda
      card.N.hat <- length(N.hat)
      
      slice.p <- (P.full.frob > threshold) # length p, true or false
      slice.pM <- rep(slice.p, each=M)
      
      A.X.eff <- A.X[, slice.pM] # select only the columns of A.X with selected features
      
      # Step 3. B check, related to lambda and epsilon
      B.check <- P.full
      for(k in 1:(p-1)){
        if(P.full.frob[k] <= threshold){
          # block not large enough: filter out this block
          B.check[(k-1)*M + (1:M), ] <- 0
        }
      }
      
      SCV.single <- rep(NA, K)
      
      for(k in 1:K){
        if(card.N.hat > 0){ # Circumstance A. if |N.hat| >=1
          
          # Step A4: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n.2/K) + 1, floor(k*n.2/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.X.cv <- A.X.eff[fold.k.ind, ]
          A.X.train <- A.X.eff[-fold.k.ind, ]
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step A5: calculate B.tilde from a ridge regression
          B.tilde <- solve(t(A.X.train) %*% A.X.train + 0.1*diag(card.N.hat * M)) %*% t(A.X.train) %*% A.Y.train
          
          # Step A6: evaluate SCV for this single k
          residual <- A.Y.cv - A.X.cv %*% B.tilde # residual is fold.k.size*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat
          
        }else{ # Circumstance B. if cardinality of N.hat is 0
          # Step B4: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n.2/K) + 1, floor(k*n.2/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step B6: evaluate SCV for this single k
          residual <- A.Y.cv # residual is fold.k.size*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          SCV.single[k] <- norm(residual,"F")^2
        }
      }# end of the k'th fold
      
      # Step 7: averaging SCV over K folds for this specific lambda and epsilon
      SCV.mat[l, ind.t] <- mean(SCV.single)
      
    } # end of all thresholds (epsilon) 1:len.t
    
  } # end of lambda 1:L
  
  # Step 8: finding the minimal SCV and its corresponding lambda/epsilon
  scv.min <- which(SCV.mat == min(SCV.mat), arr.ind=T)
  index.optimal <- scv.min[dim(scv.min)[1], ]
  l.optimal <- index.optimal[1]
  ind.t.optimal <- index.optimal[2]
  lambda.optimal <- lambdas[l.optimal]
  t.optimal <- t.vec[ind.t.optimal]
  
  SCV.optimal <- SCV.mat[index.optimal]
  B.optimal <- B.full[[l.optimal]]
  
  # Process the estimated B into a neighborhood selection vector
  B.frob <- rep(0, p-1)
  for(k in 1:(p-1)){
    B.frob[k] <- norm(B.optimal[(k-1)*M + (1:M), ], "F")
  }
  
  B.j <- rep(0,p)
  for(juliet in 1:(p-1)){
    if(juliet < j){
      B.j[juliet] <- B.frob[juliet]
    }else{
      B.j[juliet+1] <- B.frob[juliet]
    }
  }
  
  B.j
}

B.control.sym <- matrix(NA, p, p)
for(j in 1:p){
  for(k in 1:p){
    B.control.sym[j,k] <- min(B.mat.control[j,k], B.mat.control[k,j]) # equivalent to AND
  }
}

thres.vec <- seq(0, max(B.control.sym)*1.01, length.out=500)

for(i in 1:500){
  threshold <- thres.vec[i]
  G.control.sym <- matrix(0, p, p)
  for(j in 1:p){
    for(k in 1:p){
      G.control.sym[j,k] <- as.numeric(B.control.sym[j,k] > threshold)
      if(j==k) G.control.sym[j,k] <- 0
    }
  }
  if(sum(G.control.sym) <= spar.ctrl * p*(p-1)) break
}

time.end.control <- proc.time()
runtime.control <- (time.end.control - time.start.control)[3]


stopCluster(cl)

####################################
#  PART 5: SAVE ADJACENCY MATRIX   #
####################################
save(G.ADHD.sym,    file = paste(save.path, "G.spa.ctrl.ADHD.RData",    sep="/"))
save(G.control.sym, file = paste(save.path, "G.spa.ctrl.control.RData", sep="/"))


####################################
#      PART 6: SAVE RUNTIME        #
####################################
# real.data.runtime <- c(runtime.ADHD, runtime.control, runtime.fpc)
# save(real.data.runtime, file=paste(runtime.path, "real.data.Runtime.Rdata",sep="/"))

time.end <- proc.time()
time.run <- (time.end - time.start)[3]
print(time.run)