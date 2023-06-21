##################################
#  Analyzing adj mat of real data
##################################

p <- 116

library(igraph)

load("Desktop/GroupLasso/Real data analysis/Results/ABIDE/G.spa.ctrl.autism.RData")
load("Desktop/GroupLasso/Real data analysis/Results/ABIDE/G.spa.ctrl.control.RData")

# Sparsity Analysis
cat("==============================\nSparsity of adjacency matrices\n==============================\n")
spar.autism.sym <- sum(G.autism.sym)/(p*(p-1))
spar.control.sym <- sum(G.control.sym)/(p*(p-1))
cat(paste("Autism Group:  ", round(spar.autism.sym, 3)*100, "%\n", sep=""))
cat(paste("Control Group: ", round(spar.control.sym, 3)*100, "%\n", sep=""))


#####
# Visualize matrix
# par(mfrow=c(2,2))
# image(t(G.autism.and[nrow(G.autism.and):1,]), axes=FALSE, main="autism Group, AND")
# image(t(G.autism.or[nrow(G.autism.or):1,]), axes=FALSE, main="autism Group, OR")
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

hub.detect(G.autism.sym, 0.05)
hub.detect(G.control.sym, 0.05)

sum(rowSums(G.autism.sym) == 0)
sum(rowSums(G.control.sym) == 0)
sum(rowSums(G.autism.sym) == 1)
sum(rowSums(G.control.sym) == 1)
sum(rowSums(G.autism.sym) == 2)
sum(rowSums(G.control.sym) == 2)
sum(rowSums(G.autism.sym) == 3)
sum(rowSums(G.control.sym) == 3)
sum(rowSums(G.autism.sym) == 4)
sum(rowSums(G.control.sym) == 4)
sum(rowSums(G.autism.sym) == 5)
sum(rowSums(G.control.sym) == 5)
sum(rowSums(G.autism.sym) >= 6)
sum(rowSums(G.control.sym) >= 6)

sd(rowSums(G.autism.sym))
sd(rowSums(G.control.sym))

max(rowSums(G.autism.sym))
max(rowSums(G.control.sym))
