###########################
## Preliminary functions ##
###########################

# (a) Break A.X into A.X.Grp[[1]] to A.X.Grp[[p]]
A.X.into.groups <- function(A.X,M){
  # Input: A.X, n*pM matrix; M, parameter
  # Output: list A.X.Grp, length p, each n*M matrix
  n <- nrow(A.X)
  p <- ncol(A.X)/M
  A.X.Grp <- list()
  for(j in 1:p){
    A.X.Grp[[j]] <- A.X[,(j-1)*M + (1:M)]
  }
  return(A.X.Grp)
}

# (b) Ruiz scaling function
ruiz.scale <- function(A.in, tol=1e-4, max.iter=200){
  # Input:
  #   A.in, a m*n matrix
  #
  # Output:
  #   A, a scaled m*n matrix, A = D * A.in * E
  #   D, left diagonal
  #   E, right diagonal
  
  # Init
  m <- nrow(A.in); n <- ncol(A.in)
  d1 <- rep(1,m); d2 <- rep(1,n)
  A <- A.in
  
  res.r.vec <- numeric(0)
  res.c.vec <- numeric(0)
  
  # scaling
  for(ruiz.iter in 1:max.iter){
    row.max <- rep(NA, m)
    col.max <- rep(NA, n)
    
    for(i in 1:m){
      row.max[i] <- norm(matrix(A[i,], nrow=1, ncol=n), type="I")
      d1[i] <- d1[i] / sqrt(row.max[i])
    }
    for(j in 1:n){
      col.max[j] <- norm(matrix(A[,j], nrow=m, ncol=1), type="I")
      d2[j] <- d2[j] / sqrt(col.max[j])
    }
    
    
    D <- diag(d1)
    E <- diag(d2)
    A <- D %*% A.in %*% E
    
    # Check stopping criterion
    res.r <- max(row.max)/min(row.max) - 1
    res.r.vec <- c(res.r.vec, res.r)
    res.c <- max(col.max)/min(col.max) - 1
    res.c.vec <- c(res.c.vec, res.c)
    if(res.r<tol & res.c<tol){
      break
    }
  }
  result <- list(A=A, D=D, E=E, d=d1, e=d2, res.r=res.r.vec, res.c=res.c.vec,
                 ruiz.iter=ruiz.iter)
  return(result)
}

# (c) Group soft thresholding
grp.soft.thres <- function(V.Grp, d, lambda, rho){
  # V.Grp: a group of p matrices (Vk), each dim M * M
  # d: vector, length pM.
  # lambda: penalty param of group lasso
  # rho: penalty param of augmented Lagrangian
  
  M <- nrow(V.Grp[[1]])
  p <- length(V.Grp)
  D.Grp <- list()
  P.Grp <- list() # output group of p matrices (Pk), each dim M * n
  
  for(k in 1:p){
    D.Grp[[k]] <- diag(d[(k-1)*M+(1:M)])
    left.plus <- matrix(0, M, M)
    norm.F <- norm(D.Grp[[k]] %*% V.Grp[[k]], "F")
    for(l in 1:M){
      left.plus[l,l] <- max(0, 1 - lambda / rho * D.Grp[[k]][l,l]^2 / norm.F )
    }
    P.Grp[[k]] <- left.plus %*% V.Grp[[k]]
  }
  
  return(P.Grp)
}

# (d) Objective function
objective.function <- function(A.X, A.Y, d, B, lambda){
  # A.X: matrix, n * pM
  # A.Y: matrix, n * M
  # d: vector, length pM
  # P.Grp: list of length p, each with matrix M * M
  # Q: matrix, pM * M
  # lambda: penalty parameter
  
  n <- nrow(A.X)
  M <- ncol(A.Y)
  p <- ncol(A.X) / M
  
  D <- diag(d)
  term.1 <- norm(A.X %*% D %*% B - A.Y, "F")^2 / (2*n)
  term.2 <- 0
  for(k in 1:p){
    B.k <- B[(k-1)*M + (1:M), ]
    Dkk <- diag(d[(k-1)*M + (1:M)])
    term.2 <- term.2 + lambda * norm(Dkk %*% B.k ,"F")
  }
  obj <- term.1 + term.2
  return(obj)
}

# (e) kappa plus: effective condition number
kappa.plus <- function(A){
  
  sv <- svd(A)$d
  n <- length(sv)
  for(i in 1:(n-1)){
    if(sv[i+1] / sv[i] < 1e-6){
      break
    }
  }
  kappa <- sv[1] / sv[i]
  return(kappa)
}

# (f) lambda sup
lambda.sup <- function(A.X, A.Y){
  n <- nrow(A.X)
  M <- ncol(A.Y)
  p <- ncol(A.X)/M
  candidates <- rep(0,p)
  for(j in 1:p){
    A.X.j <- A.X[,(j-1)*M + (1:M)]
    candidates[j] <- norm(t(A.X.j) %*% A.Y,"F") / n
  }
  return(max(candidates))
}

# (g1) find a global lambda sup given FPC score (for gX)
lambda.sup.global.gX <- function(fpc.score, p){

  n <- nrow(fpc.score)
  M <- ncol(fpc.score) / p

  A.Y <- matrix(NA, n, M)
  A.X <- matrix(NA, n, (p-1)*M)
  l.max <- rep(0,p)
  for(j in 1:p){
    jth.range <- (((j-1)*M+1) : (j*M))
    for(i in 1:n){
      A.Y[i,] <- fpc.score[i, jth.range]
      A.X[i,] <- fpc.score[i, -jth.range]
    }
    l.max[j] <- lambda.sup(A.X, A.Y)
  }
  lambda.max <- max(l.max)
  return(lambda.max)
}

# (g2)
lambda.sup.global.gY <- function(h, M){
  
  p <- dim(h)[2]
  
  l.max <- rep(0,p)
  
  for(j in 1:p){
    A.Y <- FPC.score.j(h, j, M)$Y # n * M
    A.X <- FPC.score.j(h, j, M)$X # n * pM
    l.max[j] <- lambda.sup(A.X, A.Y)
  }
  lambda.max <- max(l.max)
  return(lambda.max)
}

###########################
##     Main Function     ##
###########################

ADMM.grplasso <- function(A.X, A.Y, d, lambda,
                          rho.init=1, maxiter=2000, mu=10, tau=2,
                          P.in, Q.in, U.in, tol.rel=1e-4, tol.abs=1e-4){
  
  # Input:
  #   A.X, n*pM matrix. For convenience we denote dimension by p instead of p-1
  #   A.Y, n*M matrix
  #   d, a pM long vector
  
  # 1.1. Reading in the parameters
  n <- nrow(A.Y)
  M <- ncol(A.Y)
  p <- ncol(A.X) / M # corresponds to p-1 in theory
  
  # 1.2. Recover matrices used in computation
  D <- diag(d)
  
  AX.D <- A.X %*% D
  DXXD.n <- (t(AX.D) %*% AX.D) / n
  DXY.n <- (t(AX.D) %*% A.Y)/n
  
  # 1.3. Setting initial values for P,Q,U
  P.old <- P.in
  Q.old <- Q.in
  U.old <- U.in
  
  # 1.4. Initialize output quantities
  prim.resid.vec <- numeric(0)
  dual.resid.vec <- numeric(0)
  rho.path <- rho.init
  tol.prim.path <- numeric(0)
  tol.dual.path <- numeric(0)
  obj.func.path <- numeric(0)
  
  # 1.5. Cached quantities -- initial
  # 1.5.1. Initial Vk's
  V <- Q.old - U.old # pM * M
  V.Grp <- list()
  for(k in 1:p){
    V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * M
  }
  
  # 1.5.2. indicator of rho's change
  rho <- rho.init
  
  ##############################
  ## 2. Main part: Iterations ##
  ##############################
  for(iter in 1:maxiter){
    
    if(iter==1){
      rho.changed <- T # For iteration 1, we need it True for Q update
    }
    
    #############################
    ## 7.1. Updating Variables ##
    #############################
    # 7.1.1. Update P with group soft thresholding
    P.new <- matrix(NA, nrow=p*M, ncol=M)
    P.Grp <- grp.soft.thres(V.Grp, d, lambda, rho) # Group of p, each M * M
    for(k in 1:p){
      P.new[(k-1)*M + (1:M), ] <- P.Grp[[k]]
    }
    # 7.1.2. Update Q
    if(rho.changed){
      # if rho not changed, this term stays the same
      left.inv <- solve(DXXD.n + rho * diag(p*M))
    }
    Q.new <- left.inv %*% (DXY.n + rho * P.new + rho * U.old)
    
    # 7.1.3. Update U
    U.new <- U.old + P.new - Q.new
    
    
    #####################################
    ## 7.2. Checking Stopping Criteria ##
    #####################################
    
    # 7.2.1. Calculate primal and dual residuals
    prim.resid <- norm( (P.new - Q.new) , "F")
    dual.resid <- norm( rho * (Q.new - Q.old) , "F")
    prim.resid.vec <- c(prim.resid.vec,prim.resid)
    dual.resid.vec <- c(dual.resid.vec,dual.resid)
    
    # 7.2.2. Calculate tolerance
    tol.prim <- tol.abs * sqrt(p*M*M) + tol.rel * max(norm(P.new, "F"), 
                                                      norm(Q.new, "F"))
    tol.dual <- tol.abs * sqrt(p*M*M) + tol.rel * norm(U.new, "F")
    # sqrt(pM^2) is because the Frob norms are in R^{pM * M}.
    
    # 7.2.3. Compare residuals with tolerance
    if(prim.resid < tol.prim & dual.resid < tol.dual){
      break
    }
    
    #####################################
    ## 7.3. Prepare for next iteration ##
    #####################################
    
    # 7.3.1. Update V
    V <- Q.new - U.new # pM * M
    V.Grp <- list()
    for(k in 1:p){
      V.Grp[[k]] <- V[(k-1)*M + (1:M), ] # each Vk dim M * n
    }
    
    # 7.3.2. Update rho for next round
    if(prim.resid > mu * dual.resid){
      rho <- tau * rho
      rho.changed <- T
    }else if(dual.resid > mu * prim.resid){
      rho <- rho / tau
      rho.changed <- T
    }else{
      rho.changed <- F
    }
    rho.path <- c(rho.path, rho)
    
    # 7.3.3. old <- new, end of one iterate
    P.old <- P.new
    Q.old <- Q.new
    U.old <- U.new
    
    
    # 7.4.1. Record objective function
    obj <- objective.function(A.X=A.X, A.Y=A.Y, d=d, B=P.new, lambda=lambda)
    obj.func.path <- c(obj.func.path, obj)
  }
  
  cat(paste("ADMM converges after ", iter," iterations.\n", sep=""))
  
  result <- list(P=P.new, Q=Q.new, U=U.new,
                 prim.res=prim.resid.vec, dual.res=dual.resid.vec,
                 rho.path=rho.path,
                 tol.prim.path=tol.prim.path, tol.dual.path=tol.dual.path,
                 obj.func.path=obj.func.path)
  return(result)
}
