# Parameters
mod <- "A"
p <- 100
array <- c(101,104,106)

if(p == 50) pee <- "050"
if(p == 100) pee <- "100"
if(p == 150) pee <- "150"

save.path <- ""
source("prec.rec.R")

n.par <- length(array)
prec.and <- 0; rec.and <- 0; prec.or <- 0; rec.or <- 0

p.and <- rep(0,n.par)
r.and <- rep(0,n.par)
p.or <- rep(0,n.par)
r.or <- rep(0,n.par)

for(i in 1:n.par){
  run.ind <- array[i]
  load(paste(save.path,"/SCV.",mod,".",pee,".RunInd",run.ind,".RData", sep=""))
  G.mat <- output.result$G.mat
  G.true <- output.result$G.true
  
  p.and[i] <- prec.rec(G.true, G.mat, "AND")$prec
  r.and[i] <- prec.rec(G.true, G.mat, "AND")$rec
  p.or[i] <- prec.rec(G.true, G.mat, "OR")$prec
  r.or[i] <- prec.rec(G.true, G.mat, "OR")$rec
  
  prec.and <- mean(p.and)
  rec.and <- mean(r.and)
  prec.or <- mean(p.or)
  rec.or <- mean(r.or)
  
  sd.prec.and <- sd(p.and)
  sd.rec.and <- sd(r.and)
  sd.prec.or <- sd(p.or)
  sd.rec.or <- sd(r.or)
}
cat(paste("Model ",mod," p=",p," Precision and Recall\n",sep=""))
cat("=================================\n")
cat(paste("AND precision: ", round(prec.and,3)," (",round(sd.prec.and,3) ,")","\n", sep=""))
cat(paste("AND recall: ", round(rec.and,3)," (",round(sd.rec.and,3) ,")","\n", sep=""))
cat(paste("OR precision: ", round(prec.or,3)," (",round(sd.prec.or,3) ,")","\n", sep=""))
cat(paste("OR recall: ", round(rec.or,3)," (",round(sd.rec.or,3) ,")","\n", sep=""))