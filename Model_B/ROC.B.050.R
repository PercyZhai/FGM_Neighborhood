###############################
##   Model_B (Model 7) ROC   ##
##           p=50            ##
###############################

func.path <- "/home/percysjzhai/GroupLasso/Functions"
working.path <- "/home/percysjzhai/GroupLasso/Model_B"
save.path <- "/home/percysjzhai/GroupLasso/Model_B/Results/p50/ROC"

# For local testing, use the following paths:
#func.path <- "~/GitHub/GroupLasso/Functions"
#working.path <- "~/GitHub/GroupLasso/Model_B"
#save.path <- "~/GitHub/GroupLasso/Model_B/Results/p50/ROC"

# Packages
library(fda)
library(matrixcalc)

# Global Parameter Settings

mu <- 15 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
p <- 50
n <- 100
tau <- 100 # number of observations
L <- 5 # number of bases used to estimate a_ij hat

# Read R files
source(paste(func.path,"prep.func.R", sep="/"))
source(paste(func.path,"ADMM.R", sep="/"))
source(paste(func.path,"grplasso.iter.R", sep="/"))
source(paste(func.path,"scv.bic.R", sep="/"))
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))
source(paste(func.path,"prec.rec.R", sep="/"))
source(paste(func.path,"auc.R", sep="/"))
source(paste(func.path,"FPCA.score.R", sep="/"))
source(paste(working.path,"M7.prelim.func.R", sep="/"))

### Reads in arguments passed in from command line
### args is a vector of strings 
args <- (commandArgs(TRUE))

for(i in 1:length(args)){
  # run each element of the vector as if passed directly to the console
  # have run.ind
eval(parse(text = args[[i]]))
}

time.start <- proc.time()


####################################
#     PART 1: DATA GENERATION      #
####################################

#    Generate Random Functions     
#      and Observation h_ijk       

# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate precision matrix and real adjacency matrix
set.seed(run.ind)
Theta <- cov.mat.model7(p, mu) # p*mu by p*mu large square matrix
G.true.7 <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true.7[i,j] <- 1
  }
}

# 1. Generating delta
delta <- rmvnorm(n, sigma = solve(Theta))

# 2. Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# 3. Fourier basis function for data generation
b.mat <- obs.fourier.bases(obs.time, mu)

# 4. Observations h_ijk
h <- array(0, c(n, p, tau))
for(i in 1:n){
  for(j in 1:p){
    h[i,j,] <- b.mat %*% matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
  }
}

# Reserved part for Zapata's method
y.list <- list()
for(j in 1:p){
  y.list[[j]] <- h[,j,]
}


####################################
#     PART 2: GAIN FPC SCORE       #
####################################

# For the use of traditional control group
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


####################################
#    PART 3: ADMM GROUP LASSO      #
#        USING FPCA BASIS          #
#        AND GET ROC CURVE         #
####################################

library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)

####################################
# Comparison: Control, Qiao, Zapata
####################################

# 3.1. Get a global lambda.max

Y <- matrix(NA, n, M)
X <- matrix(NA, n, (p-1)*M)
l.max <- rep(0,p)
for(j in 1:p){
  jth.range <- (((j-1)*M+1) : (j*M))
  for(i in 1:n){
    Y[i,] <- fpc.score[i, jth.range]
    X[i,] <- fpc.score[i, -jth.range]
  }
  l.max[j] <- lambda.sup(X,Y)
}
lambda.max <- max(l.max)
lambda.min <- 0.001
L <- 100
lambdas <- seq(lambda.max, lambda.min, length.out=L)


# 3.2. Comparison: Zapata's group graphical lasso approach
# input: directly use the observed values

n.a <- 11
alphas <- seq(0,1,length.out = n.a)

# Run across all alphas and pick the TPR/FPR that is farthest away from y=x
auc.vec <- rep(0,n.a)
TPR.zapata.list <- list()
FPR.zapata.list <- list()

for(k in 1:n.a){
  alpha <- alphas[k]
  print(paste("alpha",k,"=", alpha))
  TFPR.zapata.k <- foreach(l = 1:L, .combine="rbind", .packages="fgm") %dopar% {
    G.zapata <- fgm(y=y.list, L=5, gamma=lambdas[l], alpha=alpha)$A
    TPR.zapata.l <- prec.rec(G.true.7, G.zapata, type="AND")$TPR
    FPR.zapata.l <- prec.rec(G.true.7, G.zapata, type="AND")$FPR
    TFPR.zapata.l <- c(TPR.zapata.l, FPR.zapata.l)
    TFPR.zapata.l
  }
  
  auc.res <- auc(TFPR.zapata.k[,1], TFPR.zapata.k[,2])
  auc.vec[k] <- auc.res$AUC
  TPR.zapata.list[[k]] <- auc.res$TPR
  FPR.zapata.list[[k]] <- auc.res$FPR
}

k.opt <- which.max(auc.vec)
TPR.zapata <- TPR.zapata.list[[k.opt]]
FPR.zapata <- FPR.zapata.list[[k.opt]]



# 3.3. Comparison: Qiao's FGM, FGLasso
fpc.score.cen <- scale(fpc.score, center=T, scale=F)
est.cov.pc <- (t(fpc.score.cen) %*% fpc.score.cen) / (n-1)

TFPR.qiao <- foreach(l = 1:L, .combine="rbind", .packages="matrixcalc") %dopar% {
  G.qiao <- ProxAlg_FGM(est.cov.pc, p, M, lambdas[l])$Support # p*p support
  TPR.qiao.l <- prec.rec(G.true.7, G.qiao, type="AND")$TPR
  FPR.qiao.l <- prec.rec(G.true.7, G.qiao, type="AND")$FPR
  TFPR.qiao.l <- c(TPR.qiao.l, FPR.qiao.l)
  TFPR.qiao.l
}

TPR.qiao <- TFPR.qiao[,1]
FPR.qiao <- TFPR.qiao[,2]



# 3.4. Run ADMM and gain TPR/FPR
G.list <- list()
TPR.and <- rep(0,L)
FPR.and <- rep(0,L)
TPR.or <- rep(0,L)
FPR.or <- rep(0,L)

# default warm start ABU
A.def <- list()
for(j in 1:(p-1)) A.def[[j]] <- matrix(0,M,M)
B.def <- matrix(0.1, n, M)
U.def <- matrix(0.01, n, M)

V.array <- foreach(j = 1:p, .combine="rbind") %dopar% {
  jth.range <- (((j-1)*M+1) : (j*M))
  for(i in 1:n){
    Y[i,] <- fpc.score[i, jth.range]
    X[i,] <- fpc.score[i, -jth.range]
  }
  
  A <- A.def; B <- B.def; U <- U.def
  
  V.j <- matrix(NA, L, p)
  
  for(l in 1:L){
    lambda <- lambdas[l]
    grp.lasso.result <- ADMM.grplasso(X,Y,lambda,A.init=A, B.init=B, U.init=U)
    A <- grp.lasso.result$A
    B <- grp.lasso.result$B
    U <- grp.lasso.result$U
    
    # Process the estimated A into a neighbor selection vector
    A.mean <- rep(0, p-1)
    A.nonzero <- rep(0, p-1)
    for(juliet in 1:(p-1)){
      A.mean[juliet] <- mean(A[[juliet]])
      A.nonzero[juliet] <- (A.mean[juliet]!=0)
    }
    if(sum(A.nonzero)==0){ # all A are 0
      threshold <- 0
    }else{
      threshold <- mean(abs(A.mean[which(A.nonzero!=0)])) * 0.05
    }
    
    # Neighbor recognization
    V.jl <- rep(0, p)
    for(juliet in 1:(p-1)){
      if(juliet<j){
        if(abs(A.mean[juliet]) > threshold) V.jl[juliet] <- 1
      }else{
        if(abs(A.mean[juliet]) > threshold) V.jl[juliet + 1] <- 1
      }
    }
    
    V.j[l,] <- V.jl
  }
  
  V.j
}



# From v.array to G.list[[1 to L]]
G.list <- list()
for(l in 1:L){
  G.list[[l]] <- matrix(NA, p, p)
  for(j in 1:p){
    G.list[[l]][j,] <- V.array[(j-1)*L+l,]
  }
}


# Calculate TPR and FPR
for(l in 1:L){
  TPR.and[l] <- prec.rec(G.true.7, G.list[[l]], type="AND")$TPR
  FPR.and[l] <- prec.rec(G.true.7, G.list[[l]], type="AND")$FPR
  TPR.or[l] <- prec.rec(G.true.7, G.list[[l]], type="OR")$TPR
  FPR.or[l] <- prec.rec(G.true.7, G.list[[l]], type="OR")$FPR
}


###########################
# END OF COMPARISON PART
# STARTING FPCA PART
###########################

# 3.1. Get a global lambda.max and the vector of lambdas
#########################################################

Y <- matrix(NA, n, M)
X <- matrix(NA, n, (p-1)*M)
l.max <- rep(0,p)
for(j in 1:p){
  Y <- FPC.score.j(h, j, M)$Y
  X <- FPC.score.j(h, j, M)$X
  l.max[j] <- lambda.sup(X,Y)
}
lambda.max <- max(l.max)
lambda.min <- 0.001
L <- 100 # L: length of lambdas
lambdas <- seq(lambda.max, lambda.min, length.out=L)



# 3.2. Run ADMM and gain TPR/FPR
#################################
G.list <- list()
TPR.FPC.and <- rep(0,L)
FPR.FPC.and <- rep(0,L)
TPR.FPC.or <- rep(0,L)
FPR.FPC.or <- rep(0,L)

# default warm start ABU
A.def <- list()
for(j in 1:(p-1)) A.def[[j]] <- matrix(0,M,M)
B.def <- matrix(0.1, n, M)
U.def <- matrix(0.01, n, M)

# Parallel computing for V.array
V.array <- foreach(j = 1:p, .combine="rbind", .packages="fda") %dopar% {
  
  # load values for X and Y under this specific j
  Y <- FPC.score.j(h, j, M)$Y
  X <- FPC.score.j(h, j, M)$X
  
  # init ABU
  A <- A.def; B <- B.def; U <- U.def
  # init V.j
  V.j <- matrix(NA, L, p)
  
  for(l in 1:L){
    lambda <- lambdas[l]
    grp.lasso.result <- ADMM.grplasso(X,Y,lambda,A.init=A, B.init=B, U.init=U)
    A <- grp.lasso.result$A
    B <- grp.lasso.result$B
    U <- grp.lasso.result$U
    
    # Process the estimated A into a neighborhood selection vector
    A.mean <- rep(0, p-1)
    A.nonzero <- rep(0, p-1)
    for(juliet in 1:(p-1)){
      A.mean[juliet] <- mean(A[[juliet]])
      A.nonzero[juliet] <- (A.mean[juliet]!=0)
    }
    if(sum(A.nonzero)==0){ # all A are 0
      threshold <- 0
    }else{
      threshold <- mean(abs(A.mean[which(A.nonzero!=0)])) * 0.05
    }
    
    # Neighborhood selection
    V.jl <- rep(0, p)
    for(juliet in 1:(p-1)){
      if(juliet<j){
        if(abs(A.mean[juliet]) > threshold) V.jl[juliet] <- 1
      }else{
        if(abs(A.mean[juliet]) > threshold) V.jl[juliet + 1] <- 1
      }
    }
    
    V.j[l,] <- V.jl
  }
  
  V.j
}

# Get adjacency matrix under different lambdas
# From v.array to G.list[[1 to L]]
G.list <- list()
for(l in 1:L){
  G.list[[l]] <- matrix(NA, p, p)
  for(j in 1:p){
    G.list[[l]][j,] <- V.array[(j-1)*L+l,]
  }
}


# Calculate TPR and FPR
for(l in 1:L){
  TPR.FPC.and[l] <- prec.rec(G.true.7, G.list[[l]], type="AND")$TPR
  FPR.FPC.and[l] <- prec.rec(G.true.7, G.list[[l]], type="AND")$FPR
  TPR.FPC.or[l] <- prec.rec(G.true.7, G.list[[l]], type="OR")$TPR
  FPR.FPC.or[l] <- prec.rec(G.true.7, G.list[[l]], type="OR")$FPR
}


stopCluster(cl)


####################################
#      PART 4: SAVE ROC DATA       #
####################################

roc <- cbind(TPR.and, FPR.and, TPR.or, FPR.or, TPR.qiao, FPR.qiao, TPR.zapata, FPR.zapata, TPR.FPC.and, FPR.FPC.and, TPR.FPC.or, FPR.FPC.or)
# the following code is to remove the abnormal spikes on ROC due to systematic errors
for(l in 1:6){
  len <- length(roc[,2*l-1])
  if(abs(roc[len,2*l-1])<=0.005){
    for(i in len:2){
      if(abs(roc[i, 2*l-1])<=0.005 & roc[i-1, 2*l-1]>0.005) break
    }
    roc[i:len, 2*l-1] <- 1
    roc[i:len, 2*l] <- 1
  }
}
save(roc, file=paste(save.path,"/ROC.B.050.RunInd",run.ind,".Rdata",sep=""))


#########################################################
time.end <- proc.time()
time.run <- (time.end - time.start)[3]
print(time.run)