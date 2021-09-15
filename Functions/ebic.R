# EBIC Function given data and A.hat
EBIC <- function(X,Y,beta.hat){ # beta.hat: pM2 dimensions
  
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  
  # Extract Sparsity Pattern
  sp.pattern.pM2 <- (beta.hat != 0)
  select <- rep(c(rep(T,M), rep(F,M^2-M)),p)
  sp.pattern.pM <- sp.pattern.pM2[select] # From pM^2 to pM
  select.A <- rep(c(T, rep(F,M-1)),p)
  sp.pattern.p <- sp.pattern.pM[select.A] # From pM to p
  p.train <- sum(sp.pattern.p) # a vital parameter to be used
  
  # Calculate Sigma.hat
  Z <- X.kprod(X,M)
  AX.long <- Z%*%beta.hat
  AX.n.by.M <- matrix(nrow=n, ncol=M)
  for(i in 1:n){
    for(j in 1:M){
      AX.n.by.M[i,j] <- AX.long[(i-1)*M + j]
    }
  }
  epsilon.hat <- Y - AX.n.by.M # n by M
  Sigma.hat <- emp.cov(epsilon.hat)
  Omega.hat <- solve(Sigma.hat) # inverse
  
  # First Term
  #EBIC1 <-   + 
  EBIC1 <- norm(epsilon.hat,"F")^2 / 0.5 # or mean(eig()$values)
  
  # Second Term
  EBIC2 <- p.train * log(n)
  
  # Third Term
  gamma <- max(0, 1 - log(n) / (2*log(p)) )
  EBIC3 <- 2*gamma*log(choose(p,p.train))
  
  # Final result
  return(list(EBIC=EBIC1 + EBIC2, EBIC1=EBIC1, EBIC2=EBIC2, EBIC3=EBIC3))
}



####################################
##     Tuning Parameter lambda    ##
##             Using              ##
##          Extended BIC          ##
####################################


group.lasso.EBIC.ADMM <- function(X, Y, lambdas, udf.warmstart=F, A.in, B.in, U.in){
  # Input: X and Y are both in matrix form. X: n*pM, Y: n*M
  #lambdas: a vector of lambdas to be cross validated
  # Objective: Select a best lambda using EBIC, and output its result
  
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  L <- length(lambdas)
  
  df <- rep(0,L)
  
  # Init warmstart
  if(udf.warmstart){ # if the user applies a user defined warmstart
    A <- A.in; B <- B.in; U <- U.in
  }
  else{
    A <- list()
    for(j in 1:p) A[[j]] <- matrix(0,M,M)
    B <- matrix(0.1, n, M)
    U <- matrix(0.01, n, M)
  }
  
  # Calculate Sigma.hat
  #grp.result <- group.lasso.iter(X,Y)
  #Sigma.hat <- grp.result$Sigma.hat
  #grp.beta <- grp.result$coef
  #grp.lambda <- grp.result$lambda.opt
  
  # Start searching over all lambdas
  ext.bic <- rep(0,L)
  beta.path <- numeric(0)
  ebic1 <- rep(0,L); ebic2 <- rep(0,L); ebic3 <- rep(0,L);
  A.list <- list()
  for(l in 1:L){
    
    print(paste("lambda =",lambdas[l]))
    
    # First run ADMM on the full dataset and get a sparsity pattern (group selection)
    fit <- ADMM.grplasso(X,Y,lambdas[l], A.init=A, B.init=B, U.init=U)
    beta <- fit$coefficients
    beta.path <- cbind(beta.path,beta)
    sp.pattern.Z <- (beta != 0)
    select <- rep(c(rep(T,M), rep(F,M^2-M)),p)
    sp.pattern <- sp.pattern.Z[select] # From pM^2 to pM
    select.A <- rep(c(T, rep(F,M-1)),p)
    sp.pattern.A <- sp.pattern[select.A] # From pM to p
    seq.A <- TF.to.seq(sp.pattern.A) # a vector with length p.train
    # Warm start for next round (will affect the error path. No warm start)
    A <- fit$A; B <- fit$B; U <- fit$U
    A.list[[l]] <- A
    
    # Calculating sparsity patterns
    p.train <- sum(sp.pattern.Z)/M^2 # p.cv has the same value
    df[l] <- p.train
    
    # EBIC function
    ext.bic[l] <- EBIC(X,Y,beta)$EBIC
    ebic1[l] <- EBIC(X,Y,beta)$EBIC1
    ebic2[l] <- EBIC(X,Y,beta)$EBIC2
    ebic3[l] <- EBIC(X,Y,beta)$EBIC3
    
    # ABU for warmstart in larger loop
    if(l==1){
      A.out <- A; B.out <- B; U.out <- U
      error.min <- ext.bic[1]
    }
    if(l>=2){
      if(ext.bic[l]<error.min){
        A.out <- A; B.out <- B; U.out <- U
      }
    }
  }
  
  lambda.opt <- lambdas[which.min(ext.bic)]
  df.opt <- df[which.min(ext.bic)]
  beta.hat <- beta.path[,which.min(ext.bic)]
  
  print(paste("EBIC selects lambda =",lambda.opt,"as optimal in this range, with d.f.",df.opt,"out of",p))
  #fit <- ADMM.grplasso(X,Y,lambda.opt, A.init=A, B.init=B, U.init=U)
  #beta.hat <- fit$coefficients
  return(list(lambda=lambda.opt, beta=beta.hat, ebic.path=ext.bic, df.path=df, df.opt=df.opt,
              A.opt=A.out, B.opt=B.out, U.opt=U.out, ebic1=ebic1, ebic2=ebic2, ebic3=ebic3, A.list=A.list))
}
