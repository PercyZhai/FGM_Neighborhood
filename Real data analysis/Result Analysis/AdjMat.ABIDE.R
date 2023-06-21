##################################
#  Analyzing adj mat of real data
##################################

p <- 116

library(igraph)

load("Desktop/GroupLasso/Real data analysis/Results/ABIDE/G.mat.autism.RData")
load("Desktop/GroupLasso/Real data analysis/Results/ABIDE/G.mat.control.RData")

G.autism.and    <- matrix(0, p, p)
G.autism.or     <- matrix(0, p, p)
G.control.and <- matrix(0, p, p)
G.control.or  <- matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    if(G.mat.autism[i,j]==1 & G.mat.autism[j,i]==1){
      G.autism.and[i,j] <- 1
      G.autism.and[j,i] <- 1
    }
    if(G.mat.autism[i,j]==1 | G.mat.autism[j,i]==1){
      G.autism.or[i,j] <- 1
      G.autism.or[j,i] <- 1
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
spar.autism.and <- sum(G.autism.and)/(p*(p-1))
spar.autism.or <- sum(G.autism.or)/(p*(p-1))
spar.control.and <- sum(G.control.and)/(p*(p-1))
spar.control.or <- sum(G.control.or)/(p*(p-1))
cat(paste("Autism Group AND:  ", round(spar.autism.and, 3)*100, "%\n", sep=""))
cat(paste("Autism Group OR :  ", round(spar.autism.or , 3)*100, "%\n", sep=""))
cat(paste("Control Group AND: ", round(spar.control.and, 3)*100, "%\n", sep=""))
cat(paste("Control Group OR : ", round(spar.control.or , 3)*100, "%\n", sep=""))


#####
# Visualize matrix
# par(mfrow=c(2,2))
# image(t(G.autism.and[nrow(G.autism.and):1,]), axes=FALSE, main="Autism Group, AND")
# image(t(G.autism.or[nrow(G.autism.or):1,]), axes=FALSE, main="Autism Group, OR")
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

hub.detect(G.autism.and, 0.03)
#hub.detect(G.autism.or, 0.05)
hub.detect(G.control.and, 0.03)
#hub.detect(G.control.or, 0.05)
