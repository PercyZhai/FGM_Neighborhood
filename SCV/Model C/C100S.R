###############################
##        Model_C  SCV       ##
##            p=100          ##
###############################
p <- 100


func.path <- ""
save.path <- ""
runtime.path <- ""

# Packages
library(fda)
library(matrixcalc)
library(MASS)
library(mvtnorm)

# Global Parameter Settings
mu <- 5 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
n <- 100
tau <- 100 # number of observations
tol.abs <- 1e-4  # Tolerance (absolute) in ADMM
tol.rel <- 1e-4  # Tolerance (relative) in ADMM
L <- 100 # number of lambdas
K <- 5 # number of folds of SCV
t.vec <- c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0) # vector of t, where threshold epsilon = t * lambda
len.t <- length(t.vec)

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"C.prelim.func.R", sep="/")) # For Model C generation
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))   # FGLasso

### Reads in arguments passed in from command line
### args is a vector of strings 
args <- (commandArgs(TRUE))

for(i in 1:length(args)){
  # run each element of the vector as if passed directly to the console
  # have run.ind
  eval(parse(text = args[[i]]))
}

total.time.start <- proc.time()


####################################
#     PART 1: DATA GENERATION      #
####################################

#    Generate Random Functions     
#      and Observation h_ijk       

# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate precision matrix and real adjacency matrix
set.seed(run.ind)
seed <- run.ind
G.true <- edge.gen(p=p, seed=seed)
G.list <- edge.sets(G.true, M=mu, seed=seed)
Omega.list <- Omega.gen(G.list, seed=seed)
Sigma <- Sigma.mod.C(Omega.list)

tr.Sigma <- rep(0,mu)
for(l in 1:mu){
  Sigma.l <- Sigma[(l-1)*p + 1:p, (l-1)*p + 1:p]
  tr.Sigma[l] <- sum(diag(Sigma.l))
}


# 1. Generating delta
delta <- rmvnorm(n, sigma = Sigma)

# 2. Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# 3. Fourier basis function
b.mat <- obs.fourier.bases(obs.time, mu)

# 4. Observations h_ijk
h <- array(0, c(n, p, tau))
for(i in 1:n){
  for(j in 1:p){
    delta.ij <- matrix(NA, mu, 1)
    for(l in 1:mu){
      delta.ij[l,1] <- delta[i, (l-1)*p+j]
    }
    h[i,j,] <- b.mat %*% delta.ij + matrix(rnorm(tau, 0, sqrt(0.05*(sum(tr.Sigma))/p)), ncol=1)
  }
}

# Reserved part for PSKL Method
y.list <- list()
for(j in 1:p){
  y.list[[j]] <- h[,j,]
}
####################################
#     PART 2: GAIN FPC SCORE       #
####################################

# For the use of gX group
time.start.fpc <- proc.time()
fpc.score <- numeric(0)
for(j in 1:p){
  obs.val.matrix <- matrix(0, nrow=tau, ncol=n)
  for (i in c(1:n)){
    obs.val.vec <- as.vector(h[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
  # Construct a functional data object from the observation
  # bspline basis is the default setting of this function
  # It does not mean that the basis function is bspline!
  fd.object.array <- Data2fd(argvals=obs.time, y=obs.val.matrix, basisobj=bspline.basis)
  # FPCA process
  fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
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
#cl <- makeCluster(28)

registerDoParallel(cl)


#######################################
#   PART 3: CACHE FOR MULTIPLE RUNS   #
#######################################

# 3.1. Get Ruiz scaling
A.Y <- matrix(NA, n, M)
A.X <- matrix(NA, n, (p-1)*M)

# d.array <- foreach(j = 1:p, .combine="rbind") %dopar% {
#   jth.range <- (((j-1)*M+1) : (j*M))
#   for(i in 1:n){
#     A.Y[i,] <- fpc.score[i, jth.range]
#     A.X[i,] <- fpc.score[i, -jth.range]
#   }
#   
#   # Ruiz 
#   d.j <- ruiz.scale(t(A.X) %*% A.X)$d
#   
#   d.j
# }

d.array <- matrix(1, nrow=p, ncol=(p-1)*M)

d <- list()
norm.adj <- rep(NA,p)
for(k in 1:p){
  d[[k]] <- d.array[k,]
  norm.adj[k] <- norm(d[[k]],"2")
}

# 3.2. Get maximum lambda

# 3.2.1. maximum lambda for gX
lambda.max.gX <- lambda.sup.global.gX(fpc.score, p)
lambdas.gX <- seq(lambda.max.gX, 0, length.out=L)

# 3.2.2. maximum lambda for gY
#lambda.max.gY <- lambda.sup.global.gY(h, M)
#lambdas.gY <- seq(lambda.max.gY, 0, length.out=L)



################################################
##                                            ##
## PART 5. RUN ADMM-gX OVER EACH j AND lambda ##
##                                            ##
################################################

time.start.gX <- proc.time()

# default warm start PQU
P.def <- matrix(0, (p-1)*M, M)
Q.def <- matrix(0.1, (p-1)*M, M)
U.def <- matrix(0.01, (p-1)*M, M)

# G.mat is the optimal adjacency matrix to save
G.mat.gX <- matrix(NA, p, p)

# G.mat is the optimal adjacency matrix to save
G.mat.gX <- foreach(j = 1:p, .combine="rbind") %dopar% {
  
  jth.range <- (((j-1)*M+1) : (j*M))
  for(i in 1:n){
    A.Y[i,] <- fpc.score[i, jth.range]
    A.X[i,] <- fpc.score[i, -jth.range]
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
    full.result <- ADMM.grplasso(A.X = A.X, A.Y = A.Y, d = d[[j]],
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
          fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.X.cv <- A.X.eff[fold.k.ind, ]
          A.X.train <- A.X.eff[-fold.k.ind, ]
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step A5: calculate B.tilde from a ridge regression
          B.tilde <- solve(t(A.X.train) %*% A.X.train + 0.1*diag(card.N.hat * M)) %*% t(A.X.train) %*% A.Y.train
          
          # Step A6: evaluate SCV for this single k
          residual <- A.Y.cv - A.X.cv %*% B.tilde # residual is n*M matrix
          emp.cov.residual <- t(residual) %*% residual / fold.k.size # empirical covariance of residual
          SCV.single[k] <- norm(residual,"F")^2 + log(fold.k.size) * card.N.hat
          
        }else{ # Circumstance B. if cardinality of N.hat is 0
          # Step B4: slicing out kth fold as CV set, the rest as training set
          fold.k.ind <- seq(floor((k-1)*n/K) + 1, floor(k*n/K)) # index set for the kth fold of SCV
          fold.k.size <- length(fold.k.ind) # cardinality of Il under the paper's notation
          A.Y.cv <- A.Y[fold.k.ind, ]
          A.Y.train <- A.Y[-fold.k.ind, ]
          
          # Step B6: evaluate SCV for this single k
          residual <- A.Y.cv # residual is n*M matrix
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
  index.optimal <- scv.min[dim(scv.min)[1], ] # there could be multiple. Take the last one
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
  
  threshold <- t.optimal * lambda.optimal
  
  N.hat.optimal <- sum(B.frob >threshold)
  
  
  # Neighbor recognization
  G.j <- rep(0, p)
  for(juliet in 1:(p-1)){
    if(juliet < j){
      if(B.frob[juliet] > threshold)
        G.j[juliet] <- 1
    }else{
      if(B.frob[juliet] > threshold)
        G.j[juliet + 1] <- 1
    }
  }
  
  G.j
}
stopCluster(cl)

time.end.gX <- proc.time()
runtime.gX <- (time.end.gX - time.start.gX)[3]

####################################
#  PART 4: SAVE ADJACENCY MATRIX   #
####################################
output.result <- list(G.mat=G.mat.gX, G.true=G.true)
save(output.result, file=paste(save.path,"/SCV.C.100.RunInd",run.ind,".Rdata",sep=""))

####################################
#      PART 5: SAVE RUNTIME        #
####################################
scv.runtime <- c(runtime.gX, runtime.fpc)
save(scv.runtime, file=paste(runtime.path,
                             "/SCVtime.C.100.Runind",run.ind,".Rdata",sep=""))

#########################################################
total.time.end <- proc.time()
total.time.run <- (total.time.end - total.time.start)[3]
print(total.time.run)
