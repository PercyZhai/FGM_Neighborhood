source("auc.R")

# Parameters
mod.alp <- "C"
p <- 100
array <- c(1:50)

if(p==50)  pee <- "050"
if(p==100) pee <- "100"
if(p==150) pee <- "150"

save.path.2 <- ""

# function: interpolation
interp <- function(TPR, FPR, x.list){
  # input: lists; output: y.list
  TPR <- c(0,TPR,1)
  FPR <- c(0,FPR,1)
  
  
  len.TPR <- length(TPR)
  len.x <- length(x.list)
  y.list <- rep(1, len.x)
  for(i in 1:(len.x-1)){
    for(j in 1:(len.TPR-1)){
      if((x.list[i] >= FPR[j]) & (x.list[i] < FPR[j+1])){
        int.left <- j
        int.right <- j+1
        break
      }
    }
    y.list[i] <- TPR[int.left] + (x.list[i]-FPR[int.left])*(TPR[int.right] - TPR[int.left])/(FPR[int.right] - FPR[int.left])
  }
  return(y.list)
}



n.par <- length(array)
TPR.and <- 0; TPR.or <- 0; TPR.qiao <- 0; TPR.zapata <- 0; TPR.FPC.and <- 0; TPR.FPC.or <- 0

y.and <- 0; y.or <- 0; y.qiao <- 0; y.zapata <- 0; y.FPC.and <- 0; y.FPC.or <- 0

auc.and <- rep(0,n.par); auc.or <- rep(0,n.par); auc.qiao <- rep(0,n.par);
auc.zapata <- rep(0,n.par); auc.FPC.and <- rep(0,n.par); auc.FPC.or <- rep(0,n.par)

if(mod.alp=="D"){
  TPR.PSKL.and <- 0; TPR.PSKL.or <- 0
  y.PSKL.and <- 0; y.PSKL.or <- 0
  auc.PSKL.and <- rep(0,n.par); auc.PSKL.or <- rep(0,n.par)
}

x.list <- c(seq(0, 0.2, by=0.002), seq(0.22, 1, by=0.02))

if(p==50) pee <- "050"
if(p==100) pee <- "100"
if(p==150) pee <- "150"

for(i in 1:n.par){
  run.ind <- array[i]
  
  load(paste(save.path.2,"/ROC.",mod.alp,".",pee,".RunInd",run.ind,".RData", sep=""))
  auc.and[i] <- auc(roc[,1], roc[,2])$AUC
  auc.or[i] <- auc(roc[,3], roc[,4])$AUC
  auc.qiao[i] <- auc(roc[,5], roc[,6])$AUC
  auc.zapata[i] <- auc(roc[,7], roc[,8])$AUC
  auc.FPC.and[i] <- auc(roc[,9], roc[,10])$AUC
  auc.FPC.or[i] <- auc(roc[,11], roc[,12])$AUC
  if(mod.alp=="D"){
    auc.PSKL.and[i] <- auc(roc[,13], roc[,14])$AUC
    auc.PSKL.or[i] <- auc(roc[,15], roc[,16])$AUC
  }
  
  y.and <- y.and + interp(roc[,1], roc[,2], x.list) / n.par
  y.or <- y.or + interp(roc[,3], roc[,4], x.list) / n.par
  y.qiao <- y.qiao + interp(roc[,5], roc[,6], x.list) / n.par
  y.zapata <- y.zapata + interp(roc[,7], roc[,8], x.list) / n.par
  y.FPC.and <- y.FPC.and + interp(roc[,9], roc[,10], x.list) / n.par
  y.FPC.or <- y.FPC.or + interp(roc[,11], roc[,12], x.list) / n.par
  if(mod.alp=="D"){
    y.PSKL.and <- y.PSKL.and + interp(roc[,13], roc[,14], x.list) / n.par
    y.PSKL.or <- y.PSKL.or + interp(roc[,15], roc[,16], x.list) / n.par
  }
}


print(auc.and)
print(auc.or)

plot(c(0,x.list,1), c(0,y.and,1), main=paste("ROC Model ",mod.alp," p=",p,sep=""), type="l", col=2,
     xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), lwd=1, lty=2)
lines(c(0,x.list,1), c(0,y.or,1), col=3, lwd=1, lty=2)
lines(c(0,x.list,1), c(0,y.qiao,1), col=4, lty=3)
lines(c(0,x.list,1), c(0,y.zapata,1), col=6, lty=3)
lines(c(0,x.list,1), c(0,y.FPC.and,1), col=2, lty=1, lwd=1.5)
lines(c(0,x.list,1), c(0,y.FPC.or,1), col=3, lty=1, lwd=1.5)
if(mod.alp=="D"){
  lines(c(0,x.list,1), c(0,y.PSKL.and,1), col=2, lty=4, lwd=2)
  lines(c(0,x.list,1), c(0,y.PSKL.or,1), col=3, lty=4, lwd=2)
  legend("bottomright", c("AND-gX","OR-gX","FGLasso", "PSKL","AND-gY","OR-gY", "AND-PSKL basis", "OR-PSKL basis"),
        lty=c(2,2,3,3,1,1,3,3), col=c(2,3,4,6,2,3,2,3))
}else{
  legend("bottomright", c("AND-gX","OR-gX","FGLasso", "PSKL","AND-gY","OR-gY"),
        lty=c(2,2,3,3,1,1), col=c(2,3,4,6,2,3))
}
# Calculating AUC ROC
auc.and.m <- mean(auc.and)
auc.or.m <- mean(auc.or)
auc.qiao.m <- mean(auc.qiao)
auc.zapata.m <- mean(auc.zapata)
auc.FPC.and.m <- mean(auc.FPC.and)
auc.FPC.or.m <- mean(auc.FPC.or)
if(mod.alp=="D"){
  auc.PSKL.and.m <- mean(auc.PSKL.and)
  auc.PSKL.or.m <- mean(auc.PSKL.or)
}

auc.and.sd <- sd(auc.and)
auc.or.sd <- sd(auc.or)
auc.qiao.sd <- sd(auc.qiao)
auc.zapata.sd <- sd(auc.zapata)
auc.FPC.and.sd <- sd(auc.FPC.and)
auc.FPC.or.sd <- sd(auc.FPC.or)
if(mod.alp=="D"){
  auc.PSKL.and.sd <- sd(auc.PSKL.and)
  auc.PSKL.or.sd <- sd(auc.PSKL.or)
}

cat(paste("AUC ROC of Model ",mod.alp," p=",p,"\n",sep=""))
cat("====================================\n")
cat(paste("AND-gX: ", round(auc.and.m,3), " (", round(auc.and.sd,3), ")","\n", sep=""))
cat(paste("OR-gX: ", round(auc.or.m,3), " (", round(auc.or.sd,3), ")", "\n", sep=""))
cat(paste("FGLasso: ", round(auc.qiao.m,3), " (", round(auc.qiao.sd,3), ")","\n", sep=""))
cat(paste("PSKL: ", round(auc.zapata.m,3), " (", round(auc.zapata.sd,3), ")","\n", sep=""))
cat(paste("AND-gY: ", round(auc.FPC.and.m,3), " (", round(auc.FPC.and.sd,3), ")","\n", sep=""))
cat(paste("OR-gY: ", round(auc.FPC.or.m,3), " (", round(auc.FPC.or.sd,3), ")", "\n", sep=""))
if(mod.alp=="D"){
  cat(paste("AND-PSKL: ", round(auc.PSKL.and.m,3), " (", round(auc.PSKL.and.sd,3), ")","\n", sep=""))
  cat(paste("OR-PSKL: ", round(auc.PSKL.or.m,3), " (", round(auc.PSKL.or.sd,3), ")", "\n", sep=""))
}
