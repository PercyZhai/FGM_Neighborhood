# a function to calculate TP, FP, TN, FN
# and return precision, recall, TPR and FPR.
# Input:
#   G.true, the true p*p adjacency matrix
#   G.mat, the estimated p*p adjacency matrix
#   type,
#     AND: when two nodes both recognize each other as neighbors
#     OR: when either of them recognize each other as neighbors
# Output:
#   prec, precision
#   rec, recall
#   TPR and FPR, as their names

prec.rec <- function(G.true, G.mat, type=c("AND","OR")){
  
  p <- nrow(G.true)
  TP <- 0; TN <- 0; FP <- 0; FN <- 0
  
  if(type=="AND"){
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==1 & G.mat[j,i]==1)) TP <- TP + 1
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==0 & G.mat[j,i]==0)) FN <- FN + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==1 & G.mat[j,i]==1)) FP <- FP + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==0 & G.mat[j,i]==0)) TN <- TN + 1
      }
    }
  }else if(type=="OR"){
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==1 | G.mat[j,i]==1)) TP <- TP + 1
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==0 | G.mat[j,i]==0)) FN <- FN + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==1 | G.mat[j,i]==1)) FP <- FP + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==0 | G.mat[j,i]==0)) TN <- TN + 1
      }
    }
  }
  prec <- TP / (TP + FP)
  if(TP+FP==0) prec <- 1
  
  rec <- TP / (TP + FN)
  if(TP+FN==0) rec <- 0
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR))
}
