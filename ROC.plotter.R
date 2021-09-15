source("~/GitHub/GroupLasso/Functions/auc.R")

mod.list <- c("A","B","C","D")
p.list <- c(50,100,150)

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

# START

layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13,13), ncol=3, byrow=T), heights=c(4,4,4,4,3))

for(mod.alp in mod.list){
  for(p in p.list){
    
    par(mar=rep(1.5,4))
    #par(mar=c(1,1,1,1))
    
    if(mod.alp=="A") mod <- 6
    if(mod.alp=="B") mod <- 7
    if(mod.alp=="C") mod <- 5
    if(mod.alp=="D") mod <- 3
    
    if(mod.alp=="C" & p==150) array <- c(31:60)
    else array <- c(1:30)
    
    working.path.2 <- paste("~/GitHub/GroupLasso/Model_",mod.alp,sep="")
    save.path.2 <- paste("~/GitHub/GroupLasso/Model_",mod.alp,"/Results/p",p,"/ROC", sep="")
    
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
    
    if(mod.alp=="A"){
      plot(c(0,x.list,1), c(0,y.and,1), type="l", col=2,
           xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), lwd=1.5, lty=2,
           main=paste("p=", p, sep=""), cex.axis=0.5)
    }else{
      plot(c(0,x.list,1), c(0,y.and,1), type="l", col=2,
           xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), lwd=1.5, lty=2, cex.axis=0.5)
    }
    
    
    lines(c(0,x.list,1), c(0,y.or,1), col=3, lwd=1.5, lty=2)
    lines(c(0,x.list,1), c(0,y.qiao,1), col=4, lty=3, lwd=1.5)
    lines(c(0,x.list,1), c(0,y.zapata,1), col=6, lty=3, lwd=1.5)
    lines(c(0,x.list,1), c(0,y.FPC.and,1), col=2, lty=1, lwd=1.5)
    lines(c(0,x.list,1), c(0,y.FPC.or,1), col=3, lty=1, lwd=1.5)
    if(mod.alp=="D"){
      lines(c(0,x.list,1), c(0,y.PSKL.and,1), col=5, lty=4, lwd=2)
      lines(c(0,x.list,1), c(0,y.PSKL.or,1), col=8, lty=4, lwd=2)
    }
    
  }
}

par(mai=rep(0.2,4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c(2,3,4,6,2,3,5,8)
plot_type <- c(2,2,3,3,1,1,4,4)
plot_width <- c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2, 2)
legend(x = "top",inset = 0,
       legend = c("FPCA-gX, AND", "FPCA-gX, OR", "FGLasso","PSKL", "FPCA-gY, AND", "FPCA-gY, OR", "PSKL Basis, AND", "PSKL Basis, OR"), 
       col=plot_colors, lwd=plot_width, lty=plot_type, cex=1, ncol=4)

