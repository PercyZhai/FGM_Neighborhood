# Parameters
mod.alp <- "D"
p <- 150
array <- c(1:50)
n <- length(array)

if(p==50) pee <- "050"
if(p==100) pee <- "100"
if(p==150) pee <- "150"

save.path <- ""

runtime.gX <- rep(0,n)
runtime.gY <- rep(0,n)
runtime.fgl <- rep(0,n)
runtime.pskl <- rep(0,n)
runtime.pb <- rep(0,n)

for(runind in array){
  load(paste(save.path,"/ROCtime.", mod.alp,".",pee,".Runind",runind,".RData",sep=""))
  runtime.gX[runind] <- roc.runtime[1]
  runtime.gY[runind] <- roc.runtime[4]
  runtime.fgl[runind] <- roc.runtime[2]
  runtime.pskl[runind] <- roc.runtime[3]
  if(mod.alp == "D"){
    runtime.pb[runind] <- roc.runtime[6]
  }
}

cat(paste("ROC Runtime analysis, Model ",mod.alp," p=",p,"\n", sep=""))
cat("===========================================\n")
cat(paste("gX:      ", round(mean(runtime.gX),1)," (",round(sd(runtime.gX),1), ")\n", sep=""))
cat(paste("gY:      ", round(mean(runtime.gY),1)," (",round(sd(runtime.gY),1), ")\n", sep=""))
cat(paste("FGLasso: ", round(mean(runtime.fgl),1)," (",round(sd(runtime.fgl),1), ")\n", sep=""))
cat(paste("PSKL:    ", round(mean(runtime.pskl),1)," (",round(sd(runtime.pskl),1), ")\n", sep=""))
if(mod.alp == "D"){
  cat(paste("PSKL bas:", round(mean(runtime.pb),1)," (",round(sd(runtime.pb),1), ")\n", sep=""))
}