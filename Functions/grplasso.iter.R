# Y.vec and Z will be used in GroupLasso
index <- c(NA, rep(1:p, each=M^2))

# Calculate Sigma hat (matrix y or eps n*M), return M*M
emp.cov <- function(y){
  n <- nrow(y)
  res <- t(y) %*% y / n
  return(res)
}

# Find max eigenvalue
library(RSpectra)
max.eig <- function(Sigma){
  return(abs(eigs_sym(Sigma, 1, opts=list(retvec=F))$values))
}

# Extract XGj and get maximum 2-norm over j
max.2.norm <- function(X,M){
  n <- nrow(X)
  p <- ncol(X)/M
  norm.vec <- rep(0,p)
  for(j in 1:p){
    X.Gj <- X[,(j-1)*M + (1:M)]
    norm.vec[j] <- norm(X.Gj/sqrt(n), type="2")
  }
  return(max(norm.vec))
}

# Iteration
group.lasso.iter <- function(X,Y, maxiter=100, delta=0.05, lambda.accuracy=1e-3, C.lambda=1){
  # X and Y are both in matrix form. X: n*pM, Y: n*M
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  
  Y.vec <- Y.vectorize(Y)
  epsilon <- Y.vec
  epsilon.mat <- Y
  
  Z <- X.kprod(X,M)
  Z.expand <- cbind(1,Z)
  lambda.path <- numeric(0)
  beta.iter.path <- numeric(0)
  eta.path <- numeric(0)
  grp.index <- c(NA, rep(1:p, each=M^2))
  
  for(iter in 1:maxiter){
    Sigma.hat <- emp.cov(epsilon.mat)
    eta <- max.eig(Sigma.hat)
    eta.path <- c(eta.path, eta)
    lambda.iter <- 2 * max.2.norm(X,M) * sqrt(eta) * (M+2*sqrt(log(2*p/delta)))/sqrt(n) * C.lambda
    lambda.path <- c(lambda.path,lambda.iter)
    fit <- grplasso(x=Z.expand, y=Y.vec, index=grp.index, lambda=lambda.iter, model = LinReg())
    beta.iter <- fit$coefficients
    beta.iter.path <- cbind(beta.iter.path, beta.iter)
    epsilon <- Y.vec - Z.expand %*% beta.iter
    epsilon.mat <- matrix(epsilon, nrow=n, ncol=M, byrow=T)
  
    
    if(iter>=2)
      if(abs(lambda.iter - lambda.path[iter-1]) <= lambda.accuracy){
        return(list(lambda.opt=lambda.iter, lambda.path=lambda.path,
                    coef=beta.iter, coef.path=beta.iter.path, eta.path=eta.path, Sigma.hat=Sigma.hat))
      }
  }
  print("Maximum times of iteration reached.")
  return(list(lambda.opt=lambda.iter, lambda.path=lambda.path,
              coef=beta.iter, coef.path=beta.iter.path, eta.path=eta.path, Sigma.hat=Sigma.hat))
}


