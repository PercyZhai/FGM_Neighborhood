# Given an array of TPR and FPR, calculate its AUC ROC.

auc <- function(TPR, FPR){
  len <- length(TPR)
  
  # Step 1.1: Delete the excess 1's 
  if(TPR[1]==1){
    for(i in 1:(len-1)){
      if(TPR[i]==1 & TPR[i+1]!=1) break
    }
    TPR[1:i] <- 0
    FPR[1:i] <- 0
  }
  
  # Step 1.2: Delete excess 0's
  if(TPR[len]==0){
    for(i in len:2){
      if(TPR[i]==0 & TPR[i-1]!=0) break
    }
    TPR[i:len] <- 1
    FPR[i:len] <- 1
  }
  
  # Step 2.1: prepare to calculate AUC
  TPR <- c(0,TPR,1)
  FPR <- c(0,FPR,1)
  len <- length(TPR)
  
  # Step 2.2: calculate AUC ROC
  area <- 0
  for(i in 1:(len-1)){
    delta.x <- FPR[i+1]-FPR[i]
    y.lower <- TPR[i]
    y.upper <- TPR[i+1]
    area <- area + (y.upper + y.lower) * delta.x / 2
  }
  return(list(AUC=as.numeric(area), TPR=TPR[2:(len-1)], FPR=FPR[2:(len-1)]))
}