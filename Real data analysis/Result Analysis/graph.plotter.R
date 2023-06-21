file.path <- ""
graph.path <- paste(file.path, "Graph", sep="")
library(igraph)
source(paste(file.path,"AAL.info.R",sep=""))

layout <- function(list, option=c("xy","xz","yz")){
  p <- length(list)
  result <- matrix(NA, nrow=p, ncol=2)
  coord.vec <- unlist(coordinate.list) # length is 3p
  x.vec <- rep(NA, p)
  y.vec <- rep(NA, p)
  z.vec <- rep(NA, p)
  
  for(j in 1:p){
    x.vec[j] <- coord.vec[3*(j-1) + 1]
    y.vec[j] <- coord.vec[3*(j-1) + 2]
    z.vec[j] <- coord.vec[3*(j-1) + 3]
  }
  
  if(option=="xy"){
    result[,1] <- x.vec
    result[,2] <- y.vec
  }
  if(option=="xz"){
    result[,1] <- x.vec
    result[,2] <- z.vec
  }
  if(option=="yz"){
    result[,1] <- y.vec
    result[,2] <- z.vec
  }
  
  return(result)
}

graph.plotter <- function(adj.mat, coordinate.list, name.tag.vec, 
                          option=c("xy","xz","yz"), file.path, title, color){
  p <- ncol(adj.mat)
  
  # compute number of neighbors of each node
  nbr.count <- colSums(adj.mat)
  
  net <- graph_from_adjacency_matrix(adj.mat, mode="undirected") %>%
    set_vertex_attr("label", value = name.tag.vec)
  
  plot(net, vertex.label=name.tag.vec, vertex.shape = "circle",
       vertex.color = "red", vertex.frame.color = "red",
       vertex.label.font = 1, vertex.label.family = "Helvetica",
       vertex.label.cex = 0.3, vertex.label.dist = 0.4, vertex.label.degree = pi/4,
       edge.color=color,  edge.width = 1.2, vertex.size = 0.5*(nbr.count+1),
       layout=layout(coordinate.list, "yz"), edge.curved = 0.8)
  return(0)
}

save.path <- paste(file.path, "Graphs", sep="/")

source(paste(file.path, "AdjMat.ABIDE.R", sep=""))
colnames(G.autism.and) <- name.tag.vec
row.names(G.autism.and) <- name.tag.vec
graph.plotter(G.autism.and, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.SCV.Autism", "black")
#dev.off()
graph.plotter(G.control.and, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.SCV.Control", "black")
#dev.off()

source(paste(file.path, "AdjMat.ADHD.R", sep=""))
graph.plotter(G.ADHD.or, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.SCV.ADHD", "black")
# dev.off()
graph.plotter(G.control.or, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.SCV.Control", "black")
# dev.off()
# 
source(paste(file.path, "AdjMat.ABIDE.spar.R", sep=""))
graph.plotter(G.autism.sym, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.two.pct.Autism", "darkblue")
# dev.off()
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", save.path, "ABIDE.two.pct.Control", "darkblue")
# dev.off()
# 
source(paste(file.path, "AdjMat.ADHD.spar.R", sep=""))
graph.plotter(G.ADHD.sym, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.two.pct.ADHD", "darkblue")
# dev.off()
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", save.path, "ADHD.two.pct.Control", "darkblue")
# dev.off()

source(paste(file.path, "AdjMat.ADHD.spar.FGLasso.R", sep=""))
graph.plotter(G.ADHD.sym, coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.ADHD", "darkblue")
# dev.off()
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", save.path, "FGLasso.ADHD.two.pct.Control", "darkblue")
# dev.off()
