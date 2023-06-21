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
L <- 30 # number of lambdas
K <- 5 # number of folds of SCV
t.vec <- c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3, 5, 8, 10, 15, 20) # vector of t, where threshold epsilon = t * lambda
len.t <- length(t.vec)
spar.ctrl <- 0.02

time.start <- proc.time()


G.sparsity <- function(G){
  p <- nrow(G)
  diag(G) <- 0
  spar <- sum(G)/(p*(p-1))
  return(spar)
}

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

# 2.2. Get maximum lambda
# lambda.max.1 <- lambda.sup.global.gX(fpc.score.1, p)
# lambdas.1 <- seq(lambda.max.1, 0, length.out=L)
lambdas.1 <- seq(0.0624512843780588, 0.0596125896336016, length.out=L)

# lambda.max.2 <- lambda.sup.global.gX(fpc.score.2, p)
# lambdas.2 <- seq(lambda.max.2, 0, length.out=L)
lambdas.2 <- seq(0.0441836551730722, 0.0403415982015007, length.out=L)

# 2.3. Obtain estimated covariance
fpc.score.cen.1 <- scale(fpc.score.1, center=T, scale=F)
est.cov.pc.1 <- (t(fpc.score.cen.1) %*% fpc.score.cen.1) / (n.1-1)

fpc.score.cen.2 <- scale(fpc.score.2, center=T, scale=F)
est.cov.pc.2 <- (t(fpc.score.cen.2) %*% fpc.score.cen.2) / (n.2-1)


####################################
#    PART 3: ADHD Group FGLasso    #
####################################

for(lambda in lambdas.1){
  cat(paste("ADHD Group, Lambda =", lambda,"\n"))
  G.ADHD.sym <- ProxAlg_FGM(est.cov.pc.1, p, M, lambda)$Support
  spar <- G.sparsity(G.ADHD.sym)
  cat(paste("Sparsity:", spar,"\n"))
  if(spar >= spar.ctrl) break
}

####################################
#   PART 4: Control Group FGLasso  #
####################################

for(lambda in lambdas.2){
  cat(paste("Ctrl Group, Lambda =", lambda,"\n"))
  G.control.sym <- ProxAlg_FGM(est.cov.pc.2, p, M, lambda)$Support
  spar <- G.sparsity(G.control.sym)
  cat(paste("Sparsity:", spar, "\n"))
  if(spar >= spar.ctrl) break
}


####################################
#  PART 5: SAVE ADJACENCY MATRIX   #
####################################
save(G.ADHD.sym,    file = paste(save.path, "G.spa.ctrl.ADHD.FGLasso.RData",    sep="/"))
save(G.control.sym, file = paste(save.path, "G.spa.ctrl.control.FGLasso.RData", sep="/"))


time.end <- proc.time()
time.run <- (time.end - time.start)[3]
print(time.run)