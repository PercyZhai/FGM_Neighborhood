###############################
##        Model_A  Exp1      ##
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

mu <- 15 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
n <- 100
tau <- 100 # number of observations
tol.abs <- 1e-4  # Tolerance (absolute) in ADMM
tol.rel <- 1e-4  # Tolerance (relative) in ADMM
L <- 100 # number of lambdas

thres.ctrl.list <- c(5, 2, 1, 0.6, 0.3, 0.1, 0.05, 0)
n.thres <- length(thres.ctrl.list)

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"A.prelim.func.R", sep="/")) # For Model A generation
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
Theta <- cov.mat.model.A(p, mu) # p*mu by p*mu large square matrix
G.true <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true[i,j] <- 1
  }
}

# 1. Generating delta
delta <- rmvnorm(n, sigma = solve(Theta))

# 2. Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# 3. Fourier basis function for data generation
b.mat.list <- list()
for(j in 1:p){
  b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
}

# 4. Observations h_ijk
h <- array(0, c(n, p, tau))
for(i in 1:n){
  for(j in 1:p){
    h[i,j,] <- b.mat.list[[j]] %*% 
      matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
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
# cl <- makeCluster(cores, outfile="")

# Use the following line for faster direct computation
cl <- makeCluster(28)

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


################################################
##                                            ##
## PART 4. RUN ADMM-gX OVER EACH j AND lambda ##
##                                            ##
################################################

time.start.gX <- proc.time()

G.list.gX <- list()
TPR.and.gX <- rep(0, L*n.thres)
FPR.and.gX <- rep(0, L*n.thres)
TPR.or.gX <- rep(0, L*n.thres)
FPR.or.gX <- rep(0, L*n.thres)

# default warm start PQU
P.def <- matrix(0, (p-1)*M, M)
Q.def <- matrix(0.1, (p-1)*M, M)
U.def <- matrix(0.01, (p-1)*M, M)


V.array <- foreach(j = 1:p, .combine="rbind") %dopar% {
  jth.range <- (((j-1)*M+1) : (j*M))
  for(i in 1:n){
    A.Y[i,] <- fpc.score[i, jth.range]
    A.X[i,] <- fpc.score[i, -jth.range]
  }
  
  P <- P.def; Q <- Q.def; U <- U.def
  
  V.j <- matrix(NA, L, p*n.thres)
  
  for(l in 1:L){
    lambda <- lambdas.gX[l]
    grp.lasso.result <- ADMM.grplasso(A.X = A.X, A.Y = A.Y, d = d[[j]],
                                      lambda=lambda, rho.init=1,
                                      P.in=P, Q.in=Q, U.in=U,
                                      tol.abs = tol.abs, tol.rel = tol.rel,
                                      maxiter = 400)
    P <- grp.lasso.result$P
    Q <- grp.lasso.result$Q
    U <- grp.lasso.result$U
    
    # Process the estimated P into a neighborhood selection vector
    P.frob <- rep(0, p-1)
    for(k in 1:(p-1)){
      P.frob[k] <- norm(P[(k-1)*M + (1:M), ], "F")
    }
    
    # Neighbor recognization
    # Each V.jl contains V's for n.thres many thresholds (length: p n.thres)
    # V.jl = c( V.jl(1), V.jl(2), ..., V.jl(n.thres) )
    # V.jl concatenated by row -> V.j (L * pn)
    # V.j concatenated by row -> V.array (pL * pn)
    
    V.jl <- rep(0, p * n.thres)
    
    for(i.thres in 1:n.thres){
      thres.ctrl <- thres.ctrl.list[i.thres]
      threshold <- lambda * thres.ctrl
      
      for(juliet in 1:(p-1)){
        if(juliet < j){
          if(P.frob[juliet] > threshold)
            V.jl[(i.thres - 1) * p + juliet] <- 1
        }else{
          if(P.frob[juliet] > threshold)
            V.jl[(i.thres - 1) * p + juliet + 1] <- 1
        }
      }
    }
    
    V.j[l,] <- V.jl
  }
  V.j
}


# Get adjacency matrix under different lambdas
# From V.array to G.list[[1 to L, L+1 to 2L, ..., (n.thres-1)L+1 to n.thres*L]]
# Eventually G.list[[1 to L]] corresponds to the first threshold, [[L+1 to 2L]] to second, etc.

G.list <- list()
for(i.thres in 1:n.thres){
  for(l in 1:L){
    G.list[[(i.thres-1)*L + l]] <- matrix(NA, p, p)
    for(j in 1:p){
      G.list[[(i.thres-1)*L + l]][j,] <- V.array[(j-1)*L+l, ((i.thres-1)*p+1) : (i.thres * p)]
    }
  }
}

# Calculate TPR and FPR
for(i.thres in 1:n.thres){
  for(l in 1:L){
    TPR.and.gX[(i.thres-1)*L + l] <- prec.rec(G.true, G.list[[(i.thres-1)*L + l]], type="AND")$TPR
    FPR.and.gX[(i.thres-1)*L + l] <- prec.rec(G.true, G.list[[(i.thres-1)*L + l]], type="AND")$FPR
    TPR.or.gX[(i.thres-1)*L + l] <-  prec.rec(G.true, G.list[[(i.thres-1)*L + l]], type="OR" )$TPR
    FPR.or.gX[(i.thres-1)*L + l] <-  prec.rec(G.true, G.list[[(i.thres-1)*L + l]], type="OR" )$FPR
  }
}


time.end.gX <- proc.time()
runtime.gX <- (time.end.gX - time.start.gX)[3]

stopCluster(cl)


####################################
#      PART 5: SAVE ROC DATA       #
####################################

roc <- cbind(TPR.and.gX, FPR.and.gX, TPR.or.gX, FPR.or.gX)
# the following code is to remove the abnormal spikes on ROC due to systematic errors
for(l in 1:2){
  len <- length(roc[,2*l-1])
  if(abs(roc[len,2*l-1])<=0.005){
    for(i in len:2){
      if(abs(roc[i, 2*l-1])<=0.005 & roc[i-1, 2*l-1]>0.005) break
    }
    roc[i:len, 2*l-1] <- 1
    roc[i:len, 2*l] <- 1
  }
}
save(roc, file=paste(save.path,"/E1.A.100.RunInd",run.ind,".Rdata",sep=""))

####################################
#      PART 6: SAVE RUNTIME        #
####################################
e1.runtime <- c(runtime.gX, runtime.fpc)
save(e1.runtime, file=paste(runtime.path,
                             "/E1time.A.100.Runind",run.ind,".Rdata",sep=""))


#########################################################
total.time.end <- proc.time()
total.time.run <- (total.time.end - total.time.start)[3]
print(total.time.run)