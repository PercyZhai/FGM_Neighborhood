SCVBIC <- function(X,Y,beta.hat, empty=F){ # beta.hat: pM2 dimensions
  
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  
  if(empty){ # beta.hat is all zero
    epsilon.hat <- Y
    p.train <- 0
  }
  else{ # beta.hat is not empty
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
  }
  
  
  Sigma.hat <- emp.cov(epsilon.hat)
  #Omega.hat <- solve(Sigma.hat) # inverse
  
  # First Term
  term1 <- log(det(Sigma.hat)) * n
  #term1 <- log(det(Sigma.hat)) * n + sum(diag(epsilon.hat %*% solve(Sigma.hat) %*% t(epsilon.hat))) fails at high dim
  #term1 <- sum(diag(epsilon.hat %*% t(epsilon.hat)))
  
  # Second Term
  term2 <- p.train * log(n)
  
  # Final result
  return(list(SCVBIC= term1 + term2, term1=term1, term2=term2))
}

# A function converting True/False to the sequence of true entries
TF.to.seq <- function(TF.vec){
  seq.vec <- numeric(0)
  for(i in 1:length(TF.vec)){
    if(TF.vec[i]) seq.vec <- c(seq.vec,i)
  }
  return(seq.vec)
}

# A function that returns the minimum choice of lambda that guarantees full sparsity
lambda.sup <- function(X,Y){
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  candidates <- rep(0,p)
  for(j in 1:p){
    X.Gj <- X[,(j-1)*M + (1:M)]
    candidates[j] <- norm(t(X.Gj)%*%Y,"F")/n
  }
  return(max(candidates))
}

####################################
##     Tuning Parameter lambda    ##
##             Using              ##
##   Selective Cross Validation   ##
####################################


group.lasso.SCV.ADMM.oneloop <- function(X, Y, lambdas, K=5, udf.warmstart=F, A.in, B.in, U.in){
  # Input: X and Y are both in matrix form. X: n*pM, Y: n*M
          #lambdas: a vector of lambdas to be cross validated
          #K: K-fold SCV
  # Objective: Select a best lambda using SCV-BIC, and output its result
  
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  L <- length(lambdas)
  
  if(n %% K != 0) stop("K value error: n should be integer times of K")
  
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
  
  # Init warmstart for SCV
  p.train.old <- p+1
  A.cv.old <- list(); B.cv.old <- list(); U.cv.old <- list()
  
  # Start searching over all lambdas
  error <- rep(0,L)
  beta.path <- numeric(0)
  A.list <- list()
  for(l in 1:L){
    
    print(paste("lambda =",lambdas[l]))
    
    # First run ADMM on the full dataset and get a sparsity pattern (group selection)
    A.ein <- A; B.ein <- B; U.ein <- U
    
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
    
    # Cross-Validation with selected groups
    X.s <- X[,sp.pattern]
    scv.bic <- 0
    
    # Initial ABU for cv fits: warm start as well
    n.train <- n - n/K
    p.train <- sum(sp.pattern.Z)/M^2 # p.cv has the same value
    p.train.old <- 0
    
    if(p.train!=0){
      df[l] <- p.train
      
      if(p.train!= p.train.old) print("Sparsity changed!")
      
      for(k in 1:K){
        print(paste("Selected Cross Validation fold:",k,"/",K))
        
        # Dataset Separation
        index.cv <- ((n*(k-1)/K+1) : (n*k/K))
        X.cv <- X.s[index.cv,]
        X.train <- X.s[-index.cv,]
        Y.cv <- Y[index.cv,]
        Y.train <- Y[-index.cv,]
        Z.cv <- X.kprod(X.cv,M)
        
        # Warm Start
        if(p.train != p.train.old){ # if sparsity changed
          A.cv <- list()
          for(j in 1:(p.train)){
            A.cv[[j]] <- A[[seq.A[j]]]
          } # warm start for A: partition from large A
          if(p.train.old==0){ # if A.cv has yet to be non-zero
            B.cv <- B[-index.cv,] * K/(K-1)
            U.cv <- U[-index.cv,] * K/(K-1)
          }else{
            B.cv <- B.cv.old[[k]]
            U.cv <- U.cv.old[[k]]
          }
        }
        else{ # if sparsity not changed
          A.cv <- A.cv.old[[k]]
          B.cv <- B.cv.old[[k]]
          U.cv <- U.cv.old[[k]]
        }
        
        # Fit cv model and predict
        fit.cv <- ADMM.grplasso(X.train,Y.train,lambdas[l], A.init=A.cv, B.init=B.cv, U.init=U.cv)
        beta.cv <- fit.cv$coefficients
        pred.cv <- Z.cv %*% beta.cv
        
        # Update warm start for next round, if sparsity not changed
        A.cv <- fit.cv$A; B.cv <- fit.cv$B; U.cv <- fit.cv$U
        A.cv.old[[k]] <- A.cv; B.cv.old[[k]] <- B.cv; U.cv.old[[k]] <- U.cv
        
        # SCV-BIC Criterion. MSE + BIC term
        # Choose either of the two rows below. Edit row 156 as well
        # THIS scv.bic <- scv.bic + (norm(pred.cv-Y.vectorize(Y.cv) ,"2")^2 + p.train/p * M^2 * log(n.train))
        #scv.bic <- scv.bic + (norm(pred.cv-Y.vectorize(Y.cv) ,"2")^2 + log(n.train)*M^2*p.train)
        # SCV.AIC:
        #scv.bic <- scv.bic + norm(pred.cv-Y.vectorize(Y.cv) ,"2")^2 + p.train
        # SCV.BIC:
        #scv.bic <- scv.bic + norm(pred.cv-Y.vectorize(Y.cv) ,"2")^2 + p.train * log(n.train)
        # SCV.BIC:
        scv.bic <- scv.bic + SCVBIC(X.cv, Y.cv, beta.cv, empty=F)$SCVBIC
        
      }
      p.train.old <- p.train
      error[l] <- scv.bic
      if(l==1){
        A.out <- A; B.out <- B; U.out <- U
        error.min <- error[1]
      }
      if(l>=2){
        if(error[l]<error.min){
          error.min <- error[l]
          A.out <- A; B.out <- B; U.out <- U
        }
      }
    }
    else{ # if beta hat is all 0
      for(k in 1:K){
        index.cv <- ((n*(k-1)/K+1) : (n*k/K))
        X.cv <- X.s[index.cv,]
        X.train <- X.s[-index.cv,]
        Y.cv <- Y[index.cv,]
        Y.train <- Y[-index.cv,]
        Z.cv <- X.kprod(X.cv,M)
        # Choose either in the two rows below. Edit row 136 as well
        #scv.bic <- scv.bic + (norm(Y.vectorize(Y.cv) ,"2"))^2
        #scv.bic <- scv.bic + (norm(Y.vectorize(Y.cv) ,"2"))^2
        scv.bic <- scv.bic + SCVBIC(X.cv, Y.cv, 0, empty=T)$SCVBIC
      }
      error[l] <- scv.bic
      if(l==1){
        A.out <- A; B.out <- B; U.out <- U
        error.min <- error[1]
      }
      if(l>=2){
        if(error[l]<error.min){
          error.min <- error[l]
          A.out <- A; B.out <- B; U.out <- U
        }
      }
    }
  }
  
  lambda.opt <- lambdas[which.min(error)]
  df.opt <- df[which.min(error)]
  beta.hat <- beta.path[,which.min(error)]
  
  print(paste("SCV selects lambda =",lambda.opt,"as optimal in this range, with d.f.",df.opt,"out of",p))
  #fit <- ADMM.grplasso(X,Y,lambda.opt, A.init=A, B.init=B, U.init=U)
  #beta.hat <- fit$coefficients
  return(list(lambda=lambda.opt, beta=beta.hat, error.path=error, df.path=df, df.opt=df.opt,
              A.opt=A.out, B.opt=B.out, U.opt=U.out, A.list=A.list))
}


###########################
###   LAMBDA  SEARCH    ###
###    based on KKT     ###
###                     ###
###########################
#    Run multiple loops   #

group.lasso.SCV.ADMM.multiloop <- function(X, Y, lambda.min=1, lambda.accuracy=1.5, K=5, lambdas.terms=50){
  n <- nrow(X)
  M <- ncol(Y)
  p <- ncol(X)/M
  
  lambda.max <- lambda.sup(X,Y) + 1
  round <- 1
  while(lambda.max - lambda.min > lambda.accuracy){
    cat(paste("**********************\nlambda range:", lambda.max,"to",lambda.min,"\n"))
    lambdas <- exp(seq(log(lambda.max),log(lambda.min),length.out=lambdas.terms))
    if(round==1) result <- group.lasso.SCV.ADMM.oneloop(X,Y,lambdas,K)
    else result <- group.lasso.SCV.ADMM.oneloop(X,Y,lambdas,K,udf.warmstart = T, A.in=A, B.in=B, U.in=U)
    
    plot(lambdas, result$error.path, type="l", main="Error path", xlab="lambda",ylab="error")
    
    
    if((which(lambdas==result$lambda)+1)==length(lambdas)){ # The last lambda gives the smallest error
      lambda.min <- lambda.min - (lambda.max-lambda.min)/lambdas.terms
      lambda.max <- result$lambda
    }
    else{ # the optimal lambda is not lambda.min
      lambda.max <- result$lambda
      lambda.min <- lambdas[which(lambdas==result$lambda)+1]
    }
    A <- result$A.opt; B <- result$B.opt; U <- result$U.opt
    round <- round + 1
  }
  print(paste("The final choice of lambda is",lambda.max,"with d.f.", result$df.opt,"out of",p))
  return(list(beta=result$beta, lambda=lambda.max))
}

