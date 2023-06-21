##############################################
#   Function to calculate FPC score matrix   #
#   from observation 3D array h              #
#   based on FPCA basis using j'th node      #
##############################################
# To select neighborhood for each j, recalculate an FPC score matrix


FPC.score.j <- function(h, j, M=5){
  # h: 3D array h[i,j,k]
  # j: node index
  # M: number of basis
  
  n <- dim(h)[1]
  p <- dim(h)[2]
  tau <- dim(h)[3]
  obs.time <- seq(1/tau, 1, 1/tau)
  
  # Step 1. gain FPCA basis from h[,j,]
  obs.val.matrix.j <- matrix(0, nrow=tau, ncol=n)
  for (i in 1:n){
    obs.val.matrix.j[, i] <- as.vector(h[i, j, ])
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=M)
  # create functional data object
  fd.object.array.j <- Data2fd(argvals=obs.time, y=obs.val.matrix.j, basisobj=bspline.basis)
  # get fpca basis for h[,j,]
  fpca.basis.j <- pca.fd(fd.object.array.j, nharm=M)$harmonics
  # this fpca basis contains M functions
  
  # Step 2. convert h to functional data object
  # Step 3. calculate inner product on the IDENTICAL fpca.basis.j (different from traditional fpca)
  X.j <- matrix(NA, nrow=n, ncol=(p-1)*M)
  for(juliet in 1:p){
    obs.val.matrix.juliet <- matrix(0, nrow=tau, ncol=n)
    for(i in 1:n){
      obs.val.matrix.juliet[, i] <- as.vector(h[i, juliet, ])
    }
    fd.object.array.juliet <- Data2fd(argvals=obs.time, y=obs.val.matrix.juliet, basisobj=bspline.basis)
    # the fd object above contains n functions
    
    inn.prod.juliet <- inprod(fd.object.array.juliet, fpca.basis.j) # a n*M matrix
    
    if(juliet == j){
      Y.j <- inn.prod.juliet
    }else if(juliet < j){
      juliet.range <- (((juliet-1)*M+1) : (juliet*M))
      X.j[ , juliet.range] <- inn.prod.juliet
    }else{
      juliet.range <- (((juliet-2)*M+1) : ((juliet-1)*M))
      X.j[ , juliet.range] <- inn.prod.juliet
    }
  }
  
  return(list(Y=Y.j, X=X.j))
}