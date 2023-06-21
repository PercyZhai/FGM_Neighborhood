# Parameters
mod.alp <- "A"
p <- 100
array <- 1:50
n <- length(array)

if(p==50) pee <- "050"
if(p==100) pee <- "100"
if(p==150) pee <- "150"

save.path <- ""

runtime.gX <- rep(0,n)

for(runind in array){
  load(paste(save.path,"/SCVtime.", mod.alp,".",pee,".Runind",runind,".RData",sep=""))
  runtime.gX[runind] <- scv.runtime[1]
}

cat(paste("SCV Runtime analysis, Model ",mod.alp," p=",p,"\n", sep=""))
cat("===========================================\n")
cat(paste("gX:", round(mean(runtime.gX),1)," (",round(sd(runtime.gX),1), ")\n", sep=""))
