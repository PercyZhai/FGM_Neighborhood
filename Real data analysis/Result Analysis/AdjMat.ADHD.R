##################################
#  Analyzing adj mat of real data
##################################

p <- 116

library(igraph)

load("Desktop/GroupLasso/Real data analysis/Results/ADHD/G.mat.ADHD.RData")
load("Desktop/GroupLasso/Real data analysis/Results/ADHD/G.mat.control.RData")

G.ADHD.and    <- matrix(0, p, p)
G.ADHD.or     <- matrix(0, p, p)
G.control.and <- matrix(0, p, p)
G.control.or  <- matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    if(G.mat.ADHD[i,j]==1 & G.mat.ADHD[j,i]==1){
      G.ADHD.and[i,j] <- 1
      G.ADHD.and[j,i] <- 1
    }
    if(G.mat.ADHD[i,j]==1 | G.mat.ADHD[j,i]==1){
      G.ADHD.or[i,j] <- 1
      G.ADHD.or[j,i] <- 1
    }
    if(G.mat.control[i,j]==1 & G.mat.control[j,i]==1){
      G.control.and[i,j] <- 1
      G.control.and[j,i] <- 1
    }
    if(G.mat.control[i,j]==1 | G.mat.control[j,i]==1){
      G.control.or[i,j] <- 1
      G.control.or[j,i] <- 1
    }
  }
}
# up till now, the diagonals of G's are all zero.

# Sparsity Analysis
cat("==============================\nSparsity of adjacency matrices\n==============================\n")
spar.ADHD.and <- sum(G.ADHD.and)/(p*(p-1))
spar.ADHD.or <- sum(G.ADHD.or)/(p*(p-1))
spar.control.and <- sum(G.control.and)/(p*(p-1))
spar.control.or <- sum(G.control.or)/(p*(p-1))
cat(paste("ADHD Group AND:  ", round(spar.ADHD.and, 3)*100, "%\n", sep=""))
cat(paste("ADHD Group OR :  ", round(spar.ADHD.or , 3)*100, "%\n", sep=""))
cat(paste("Control Group AND: ", round(spar.control.and, 3)*100, "%\n", sep=""))
cat(paste("Control Group OR : ", round(spar.control.or , 3)*100, "%\n", sep=""))


#####
# Visualize matrix
# par(mfrow=c(2,2))
# image(t(G.ADHD.and[nrow(G.ADHD.and):1,]), axes=FALSE, main="ADHD Group, AND")
# image(t(G.ADHD.or[nrow(G.ADHD.or):1,]), axes=FALSE, main="ADHD Group, OR")
# image(t(G.control.and[nrow(G.control.and):1,]), axes=FALSE, main="Control Group, AND")
# image(t(G.control.or[nrow(G.control.or):1,]), axes=FALSE, main="Control Group, OR")

hub.detect <- function(G, thres=0.3){
  p <- nrow(G)
  hub.count <- 0
  for(i in 1:p){
    n.conn <- sum(G[i,])
    if(n.conn >= thres*p) hub.count <- hub.count + 1
  }
  return(hub.count)
}

hub.detect(G.ADHD.and, 0.03)
#hub.detect(G.ADHD.or, 0.05)
hub.detect(G.control.and, 0.03)
#hub.detect(G.control.or, 0.05)
