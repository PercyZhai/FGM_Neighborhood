##################################
#  Analyzing adj mat of real data
##################################


#load("GitHub/GroupLasso/Real data analysis/G_mat.RData")

#G.mat.1 <- G.mat[,1:176]
#G.mat.2 <- G.mat[,177:352]

load("GitHub/GroupLasso/Real data analysis/G1_2.5.RData")
load("GitHub/GroupLasso/Real data analysis/G2_2.5.RData")

G.1 <- matrix(0, 176, 176)
# AND scheme
for(i in 1:176){
  for(j in 1:176){
    if(G.mat.1[i,j]==1 | G.mat.1[j,i]==1){
      G.1[i,j] <- 1
      G.1[j,i] <- 1
    }
  }
}

G.2 <- matrix(0, 176, 176)
# AND scheme
for(i in 1:176){
  for(j in 1:176){
    if(G.mat.2[i,j]==1 | G.mat.2[j,i]==1){
      G.2[i,j] <- 1
      G.2[j,i] <- 1
    }
  }
}

# Sparsity of G1 and G2
cat("========================\n
     Sparsity of G1 and G2\n
    ========================\n")
spar.1 <- sum(G.1)/(176*176)
spar.2 <- sum(G.2)/(176*176)
print(paste("ADHD Group:    ", round(spar.1, 3)*100, "%", sep=""))
print(paste("Control Group: ", round(spar.2, 3)*100, "%", sep=""))


#####
# Visualize matrix
par(mfrow=c(1,2))
image(t(G.1[nrow(G.1):1,]), axes=FALSE, main="ADHD Group")
image(t(G.2[nrow(G.2):1,]), axes=FALSE, main="Control Group")

#image(t(G.1[nrow(G.1):1,]) - t(G.2[nrow(G.2):1,]), axes=FALSE, main="Difference")

hub.detect <- function(G, thres=0.3){
  p <- nrow(G)
  hub.count <- 0
  for(i in 1:p){
    n.conn <- sum(G[i,])-1
    if(n.conn >= thres*p) hub.count <- hub.count + 1
  }
  return(hub.count)
}

hub.detect(G.1, 0.5)
hub.detect(G.2, 0.5)
