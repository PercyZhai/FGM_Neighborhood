### Preliminary Functions ###
#############################

X.into.groups <- function(X,M){
  # Input: X, n*pM matrix; M, parameter
  # Output: list X.G, length p, each n*M matrix
  n <- nrow(X)
  p <- ncol(X)/M
  X.G <- list()
  for(j in 1:p){
    X.G[[j]] <- X[,(j-1)*M + (1:M)]
  }
  return(X.G)
}



bar.XGA <- function(X.G, A){
  # Function calculating bar{X.G * A}
  # Input: X.G, a list of length p with X.Gj;
          # A, a list with length p with A.j
  # Output: \bar{A_G X}
  p <- length(X.G)
  result <- 0
  for(j in 1:p){
    result <- result + (X.G[[j]]%*%A[[j]])/p
  }
  return(result)
}

primal.resid <- function(X.G, A, B, Y){
  # Function calculating primal residual
  # Input: X.G, a list of length p with X.Gj;
          # A, a list of length p with A.j;
          # B, a matrix
  # Output: a real number ||r||_F
  p <- length(X.G)
  # Alternative method by Boyd
  result <- sqrt(p)*norm(bar.XGA(X.G,A)-B,"F")
  return(result)
}

dual.resid <- function(rho, X.G, B.new, B.old){
  # Function calculating dual residual
  # Input: rho, the penalty parameter
          # X.G, a list of length p with X.Gj;
          # B, a matrix
  # Output: a real number ||s||_F
  p <- length(X.G)
  result.sq <- 0
  for(j in 1:p){
    s.j <- rho * t(X.G[[j]])%*%(B.new-B.old)
    result.sq <- result.sq + norm(s.j, "F")^2
  }
  return(sqrt(result.sq))
}

obj.func <- function(X,Y,beta,lambda){
  # Calculate objective function given an estimated beta under X,Y and lambda
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  Y.vec <- Y.vectorize(Y)
  Z <- X.kprod(X,M)
  Z.expand <- cbind(1,Z)
  term1 <- norm(Y.vec-Z.expand%*%beta,"2")^2/(2*n)
  term2 <- 0
  for(j in 1:p){
    term2 <- term2 + norm(beta[(2+(j-1)*M^2):(j*M^2+1)],"2")
  }
  term2 <- lambda*term2
  result <- term1 + term2
  return(result)
}


# Setting default initial values for A,B,U
#A.def <- list()
#for(j in 1:p) A.def[[j]] <- matrix(0,M,M)
#B.def <- matrix(0.1, n, M)
#U.def <- matrix(0.01, n, M)

##### Main Function #####
#########################

ADMM.grplasso <- function(X, Y, lambda, rho.init=1, maxiter=2000, mu=10, tau=2,
                          A.init=A.def, B.init=B.def, U.init=U.def, eps.rel=1e-4, eps.abs=1e-4){
  # Our main function.
  # Input: X, n*pM matrix;
          # Y, n*M matrix (underline Y in notes)
          # lambda_n
          # A.init: list of p many M*M matrices, B.init U.init: n*M matrix
  # Output: array of beta
  
  
  # Reading in the parameters
  n <- nrow(Y)
  M <- ncol(Y)
  p <- ncol(X) / M
  X.G <- X.into.groups(X, M)
  
  
  # Setting initial values for A,B,U
  A.old <- A.init
  B.old <- B.init
  U.old <- U.init
  
  # Calculating initial V
  V <- list()
  for(j in 1:p){
    V[[j]] <- X.G[[j]]%*%A.old[[j]] + B.old - bar.XGA(X.G, A.old) - U.old
  }
  
  
  # Eigendecomposition for X.Gj^T%*%X.Gj
  Q <- list()
  eig.diag <- list()
  for(j in 1:p){
    eig <- eigen(t(X.G[[j]])%*%X.G[[j]])
    Q[[j]] <- eig$vectors
    eig.diag[[j]] <- diag(eig$values)
  }
  
  # Calculating initial residuals
  rho <- rho.init
  r <- primal.resid(X.G, A.old, B.old, Y)
  s <- dual.resid(rho, X.G, matrix(0,n,M), B.old)
  
  # paths
  prim.resid.vec <- numeric(0)
  dual.resid.vec <- numeric(0)
  rho.path <- numeric(0)
  eps.prim.path <- numeric(0)
  eps.dual.path <- numeric(0)
  beta <- rep(0,p*M*M)
  beta.path <- numeric(0)
  obj.func.path <- numeric(0)
  
  # Start iteration
  for(k in 1:maxiter){
    
    # Calculate A by iterating for nu
    A.new <- list()
    #lhs <- numeric(0)
    
    for(j in 1:p){
      temp <- t(X.G[[j]]) %*% V[[j]]
      if(norm(temp, "F") <= lambda/rho)
        A.new[[j]] <- matrix(0,M,M)
      else{
        cached <- t(Q[[j]]) %*% temp
        nu.L <- 0
        nu.R <- 1e+8
        lhs.L <- 0
        lhs.R <- 1e+8
        rhs <- lambda/rho
        while(lhs.R-lhs.L>=rhs*1e-4){
          lhs.L <- nu.L * norm(solve(eig.diag[[j]]+nu.L*diag(M)) %*% cached,"F")
          lhs.R <- nu.R * norm(solve(eig.diag[[j]]+nu.R*diag(M)) %*% cached,"F")
          nu.C <- (nu.L+nu.R)/2
          lhs.C <- nu.C * norm(solve(eig.diag[[j]]+nu.C*diag(M)) %*% cached,"F")
          if(lhs.C < rhs) nu.L <- nu.C
          else if(lhs.C >= rhs) nu.R <- nu.C
        }
        nu <- nu.C
        A.new[[j]] <- solve(t(X.G[[j]])%*%X.G[[j]] + nu*diag(M)) %*% temp
      }
    }
    
    # Update XGA bar
    XGA.bar.new <- bar.XGA(X.G, A.new)
    
    # Update B
    B.new <- (Y + rho*n*(XGA.bar.new + U.old)) / (p+rho*n)
    # Update U
    U.new <- U.old + XGA.bar.new - B.new

    # Calculate V for next round
    for(j in 1:p){
      V[[j]] <- X.G[[j]]%*%A.new[[j]] + B.new - XGA.bar.new - U.new
    }
    
    # Calculate residuals
    p.res <- primal.resid(X.G, A.new, B.new, Y)
    d.res <- dual.resid(rho, X.G, B.new, B.old)
    prim.resid.vec <- c(prim.resid.vec,p.res)
    dual.resid.vec <- c(dual.resid.vec,d.res)
    
    # Update rho for next round
    if(p.res/sqrt(p)>mu*d.res/sqrt(n)) rho <- tau*rho
    else if(d.res/sqrt(n)>mu*p.res/sqrt(p)) rho<- rho/tau
    rho.path <- c(rho.path, rho)
    
    # checking stopping criteria
    s1.sq <- 0
    for(j in 1:p){
      s1.sq <- s1.sq + norm(rho * t(X.G[[j]]) %*% U.new,"F")^2
    }
    r1 <- norm(XGA.bar.new,"F")
    r2 <- norm(B.new,"F")
    s1 <- sqrt(s1.sq)
    eps.pri <- eps.abs * sqrt(p) + eps.rel * max(r1,r2)
    eps.dual <- eps.abs * sqrt(n) + eps.rel * s1
    #print(paste("Stopping criteria:", eps.pri, eps.dual))
    eps.prim.path <- c(eps.prim.path, eps.pri)
    eps.dual.path <- c(eps.dual.path, eps.dual)
    
    if(p.res<eps.pri & d.res<eps.dual){
      print(paste("ADMM stopping criteria satisfied at round",k))
      break
    }
    
    # old <- new, end of one iterate
    A.old <- A.new
    B.old <- B.new
    U.old <- U.new
    
    # Record beta path and objective function path
    for(j in 1:p)
      for(i in 1:M)
        for(k in 1:M)
          beta[(j-1)*25 + (i-1)*5 +k] <- A.new[[j]][i,k]
    beta.path <- cbind(beta.path, beta)
    obj.func.path <- c(obj.func.path, obj.func(X,Y,c(0,beta),lambda))
  }
  
  result <- list(A=A.new, B=B.new, U=U.new, p.res=prim.resid.vec, d.res=dual.resid.vec,
                 rho.path=rho.path, eps.prim.path=eps.prim.path, eps.dual.path=eps.dual.path,
                 coefficients=beta, beta.path=beta.path, obj.func.path=obj.func.path)
  return(result)
  
}