##################################
#  Analyzing adj mat of real data
##################################

p <- 116

library(igraph)

load("Desktop/GroupLasso/Real data analysis/Results/ADHD/G.spa.ctrl.ADHD.RData")
load("Desktop/GroupLasso/Real data analysis/Results/ADHD/G.spa.ctrl.control.RData")

# Sparsity Analysis
cat("==============================\nSparsity of adjacency matrices\n==============================\n")
spar.ADHD.sym <- sum(G.ADHD.sym)/(p*(p-1))
spar.control.sym <- sum(G.control.sym)/(p*(p-1))
cat(paste("ADHD Group:  ", round(spar.ADHD.sym, 3)*100, "%\n", sep=""))
cat(paste("Control Group: ", round(spar.control.sym, 3)*100, "%\n", sep=""))


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

# hub.detect(G.ADHD.and, 0.03)
# hub.detect(G.ADHD.or, 0.05)
# hub.detect(G.control.and, 0.03)
# hub.detect(G.control.or, 0.05)

edge.ADHD <- rowSums(G.ADHD.sym)
edge.control <- rowSums(G.control.sym)