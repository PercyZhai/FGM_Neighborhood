source("auc.R")

# Parameters to change
######################
mod.alp <- "D"    ####
p <- 150          ####
array <- c(1:50)  ####
n.thres <- 8      ####
######################

n <- 100
cex <- 0.8

if(p==50) pee <- "050"
if(p==100) pee <- "100"
if(p==150) pee <- "150"

save.path <- ""

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

y.FPC.and <- list(); y.FPC.or <- list()
for(i.thres in 1:n.thres){
  y.FPC.and[[i.thres]] <- rep(0, 141)
  y.FPC.or[[ i.thres]] <- rep(0, 141)
}

auc.FPC.and <- rep(0,n.par*n.thres); auc.FPC.or <- rep(0,n.par*n.thres)

if(mod.alp=="D"){
  TPR.PSKL.and <- 0; TPR.PSKL.or <- 0
  y.PSKL.and <- 0; y.PSKL.or <- 0
  auc.PSKL.and <- rep(0,n.par); auc.PSKL.or <- rep(0,n.par)
}

x.list <- c(seq(0, 0.2, by=0.002), seq(0.22, 1, by=0.02))

for(i in 1:n.par){
  run.ind <- array[i]
  
  load(paste(save.path,"/E1.",mod.alp,".",pee, ".RunInd",run.ind,".RData", sep=""))
  
  for(i.thres in 1:n.thres){
    roc.i.thres <- roc[seq(((i.thres-1)*n+1), i.thres*n) , ]
    
    auc.FPC.and[(i.thres-1)*n.par + i] <- auc(roc.i.thres[,1], roc.i.thres[,2])$AUC
    auc.FPC.or[ (i.thres-1)*n.par + i] <- auc(roc.i.thres[,3], roc.i.thres[,4])$AUC
    
    y.FPC.and[[i.thres]] <- y.FPC.and[[i.thres]] + interp(roc.i.thres[,1], roc.i.thres[,2], x.list) / n.par
    y.FPC.or[[i.thres]]  <- y.FPC.or[[ i.thres]] + interp(roc.i.thres[,3], roc.i.thres[,4], x.list) / n.par
  }
}

# print(auc.FPC.and)
# print(auc.FPC.or)

# Calculating AUC ROC
auc.FPC.and.m  <- rep(0, n.thres)
auc.FPC.or.m   <- rep(0, n.thres)
auc.FPC.and.sd <- rep(0, n.thres)
auc.FPC.or.sd  <- rep(0, n.thres)

for(i.thres in 1:n.thres){
  auc.FPC.and.m[i.thres]  <- mean(auc.FPC.and[seq(((i.thres-1)*n.par+1), i.thres*n.par)])
  auc.FPC.or.m[i.thres]   <- mean(auc.FPC.or[ seq(((i.thres-1)*n.par+1), i.thres*n.par)])
  auc.FPC.and.sd[i.thres] <- sd(  auc.FPC.and[seq(((i.thres-1)*n.par+1), i.thres*n.par)])
  auc.FPC.or.sd[i.thres]  <- sd(  auc.FPC.or[ seq(((i.thres-1)*n.par+1), i.thres*n.par)])
}

################################################################################################
################################################################################################


################################
##                            ##
##         ROC A 050 E        ##
##      AND 1468, OR 1468     ##
################################
## c(3, 2, 1, 0.6, 0.3, 0.1, 0.05, 0)
## Best: AND 0.6, OR 0.6
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1), main=paste("p=50"),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(4,8)){
#   if(i==4){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==4){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 2 ; AND",
#                         "t = 0.6 ; AND",
#                         "t = 0; AND",
#                         "t = 2 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#                         col=c(2,4,8), lwd=c(1.5, 2.5, 1.5),
#                         lty=c(rep(1,3), rep(2,3)))

# cat(paste("E1 AUC ROC, Model A, p=50\n",sep=""))
# cat("====================================\n")
# for(i in c(2,4,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,4,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }

################################
##                            ##
##         ROC A 100 E        ##
##      AND 1468, OR 1468     ##
################################
# c(5, 2, 1, 0.6, 0.3, 0.1, 0.05, 0)
# Best: AND 0.3, OR 0.6
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      main=paste("p=100"),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(5,8)){
#   if(i==5){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==4){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 2 ; AND",
#                         "t = 0.3 ; AND",
#                         "t = 0; AND",
#                         "t = 2 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,5,8,2,4,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))

# cat(paste("E1 AUC ROC, Model A, p=100\n",sep=""))
# cat("====================================\n")
# for(i in c(1,5,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(1,4,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }

################################
##                            ##
##         ROC A 150 E        ##
##      AND 1568, OR 1468     ##
################################
# thres.ctrl.list <- c(3, 2, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0.3, OR 0.6
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      main=paste("p=150"),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(5,8)){
#   if(i==5){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==4){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 2 ; AND",
#                         "t = 0.3 ; AND",
#                         "t = 0; AND",
#                         "t = 2 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,5,8,2,4,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model A, p=150\n",sep=""))
# cat("====================================\n")
# for(i in c(2,5,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,4,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC B 050 E        ##
##      AND 2358, OR 2368     ##
################################
# c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0)
# best: AND 0.05, OR 0.02
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(5,8)){
#   if(i==5){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=4, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(6,8)){
#   if(i==6){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=5, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 0.5 ; AND",
#                         "t = 0.05 ; AND",
#                         "t = 0; AND",
#                         "t = 0.5 ; OR",
#                         "t = 0.02 ; OR",
#                         "t = 0; OR"),
#        col=c(2,4,8,2,5,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))

# cat(paste("E1 AUC ROC, Model B, p=50\n",sep=""))
# cat("====================================\n")
# for(i in c(2,3,5,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,3,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC B 100 E        ##
##      AND 2368, OR 2368     ##
################################
# c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0)
# best: AND 0.02, OR 0.02
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(6,8)){
#   if(i==6){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=5, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(6,8)){
#   if(i==6){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=5, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 0.5 ; AND",
#                         "t = 0.02 ; AND",
#                         "t = 0; AND",
#                         "t = 0.5 ; OR",
#                         "t = 0.02 ; OR",
#                         "t = 0; OR"),
#        col=c(2,5,8,2,5,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model B, p=100\n",sep=""))
# cat("====================================\n")
# for(i in c(2,3,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,3,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC B 150 E        ##
##      AND 2368, OR 2368     ##
################################
# c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0)
# best: AND 0.02, OR 0.02
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(6,8)){
#   if(i==6){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=5, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(6,8)){
#   if(i==6){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=5, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 0.5 ; AND",
#                         "t = 0.02 ; AND",
#                         "t = 0; AND",
#                         "t = 0.5 ; OR",
#                         "t = 0.02 ; OR",
#                         "t = 0; OR"),
#        col=c(2,5,8,2,5,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model B, p=150\n",sep=""))
# cat("====================================\n")
# for(i in c(2,3,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,3,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC C 050 E        ##
##      AND 1368, OR 1368     ##
################################
# c(3, 1.5, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0, OR 0
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(4,8)){
#   if(i==8){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==8){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 1.5 ; AND",
#                         "t = 0.6 ; AND",
#                         "t = 0; AND",
#                         "t = 1.5 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,4,8,2,4,8), lwd=c(1.5, 1.5, 2.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model C, p=50\n",sep=""))
# cat("====================================\n")
# for(i in c(2,4,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,4,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }

################################
##                            ##
##         ROC C 100 E        ##
##      AND 2468, OR 2468     ##
################################
# c(3, 1.5, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0, OR 0
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(4,8)){
#   if(i==8){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==8){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 1.5 ; AND",
#                         "t = 0.6 ; AND",
#                         "t = 0; AND",
#                         "t = 1.5 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,4,8,2,4,8), lwd=c(1.5, 1.5, 2.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model C, p=100\n",sep=""))
# cat("====================================\n")
# for(i in c(2,4,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,4,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC C 150 E        ##
##      AND 2468, OR 2468     ##
################################
# c(3, 1.5, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0, OR 0
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[2]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[2]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(4,8)){
#   if(i==8){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==8){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 1.5 ; AND",
#                         "t = 0.6 ; AND",
#                         "t = 0; AND",
#                         "t = 1.5 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,4,8,2,4,8), lwd=c(1.5, 1.5, 2.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model C, p=150\n",sep=""))
# cat("====================================\n")
# for(i in c(2,4,6,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(2,4,6,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }

################################
##                            ##
##         ROC D 050 E        ##
##      AND 1578, OR 1478     ##
################################
# c(3, 1.5, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0.3, OR 0.6
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[1]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[1]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(5,8)){
#   if(i==5){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==4){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 3 ; AND",
#                         "t = 0.3 ; AND",
#                         "t = 0; AND",
#                         "t = 3 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,5,8,2,4,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model D, p=50\n",sep=""))
# cat("====================================\n")
# for(i in c(1,5,7,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(1,4,7,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC D 100 E        ##
##      AND 1678, OR 1478     ##
################################
# c(3, 1.5, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0.1, OR 0.6
################################
# plot(c(0,x.list,1), c(0,y.FPC.and[[1]],1),
#      type="l", col=2,
#      xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
# lines(c(0,x.list,1), c(0,y.FPC.or[[1]],1) , col=2, lty=2, lwd=1.5)
# 
# for(i in c(6,8)){
#   if(i==6){
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=3, lty=1, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
#   }
# }
# 
# for(i in c(4,8)){
#   if(i==4){
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
#   }else{
#     lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
#   }
# }
# legend("bottomright", c("t = 3 ; AND",
#                         "t = 0.1 ; AND",
#                         "t = 0; AND",
#                         "t = 3 ; OR",
#                         "t = 0.6 ; OR",
#                         "t = 0; OR"),
#        col=c(2,3,8,2,4,8), lwd=c(1.5, 2.5, 1.5),
#        lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model D, p=100\n",sep=""))
# cat("====================================\n")
# for(i in c(1,6,7,8)){
#   cat(paste("AND-gX threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(1,4,7,8)){
#   cat(paste("OR-gX  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }


################################
##                            ##
##         ROC D 150 E        ##
##      AND 1578, OR 1578     ##
################################
# c(3, 1.5, 1, 0.6, 0.3, 0.1, 0.05, 0)
# best: AND 0, OR 0.3
################################
plot(c(0,x.list,1), c(0,y.FPC.and[[1]],1),
     type="l", col=2,
     xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), lty=1, lwd=1.5, cex.axis=cex)
lines(c(0,x.list,1), c(0,y.FPC.or[[1]],1) , col=2, lty=2, lwd=1.5)

for(i in c(4,8)){
  if(i==8){
    lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=2.5)
  }else{
    lines(c(0,x.list,1), c(0,y.FPC.and[[i]],1), col=i, lty=1, lwd=1.5)
  }
}

for(i in c(5,8)){
  if(i==5){
    lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=2.5)
  }else{
    lines(c(0,x.list,1), c(0,y.FPC.or[[i]],1) , col=i, lty=2, lwd=1.5)
  }
}
legend("bottomright", c("t = 3 ; AND",
                        "t = 0.6 ; AND",
                        "t = 0; AND",
                        "t = 3 ; OR",
                        "t = 0.3 ; OR",
                        "t = 0; OR"),
       col=c(2,4,8,2,5,8), lwd=c(1.5, 1.5, 2.5, 1.5, 2.5, 1.5),
       lty=c(rep(1,3), rep(2,3)))
# 
# cat(paste("E1 AUC ROC, Model D, p=150\n",sep=""))
# cat("====================================\n")
# for(i in c(1,4,5,6,8)){
#   cat(paste("AND-gY threshold ", i,": ", round(auc.FPC.and.m[i],4),
#             " (", round(auc.FPC.and.sd[i],3), ")","\n", sep=""))
# }
# for(i in c(1,4,5,6,8)){
#   cat(paste("OR-gY  threshold ", i,": ", round(auc.FPC.or.m[i],4),
#             " (", round(auc.FPC.or.sd[i],3), ")", "\n", sep=""))
# }



# cat("====================================\n")
# cat(paste("AND-gY: ", round(auc.FPC.and.m,4), " (", round(auc.FPC.and.sd,3), ")","\n", sep=""))
# cat(paste("OR-gY: ", round(auc.FPC.or.m,4), " (", round(auc.FPC.or.sd,3), ")", "\n", sep=""))

