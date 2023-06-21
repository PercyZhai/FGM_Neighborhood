##################
#  ABIDE data reader
#  Reads in NYU data and generate n*p*T array h[i,j,t]
##################

file.path <- ""
file.list <- list.files(file.path, "NYU?")
phenotypic <- read.csv("Phenotypic_V1_0b_preprocessed1.csv")

# extract SUB_ID
sub_id <- as.numeric(sapply(strsplit(file.list,"_"),"[[",2))

n <- length(file.list)
p <- 116
tau <- 176 

h <- array(0, c(n,p,tau))
dx_group <- rep(0,n)

for(i in 1:n){
  f1 <- read.delim(paste(file.path, file.list[i], sep=""))
  f1 <- as.matrix(f1)
  # print(dim(f1)) everyone 176*116
  h[i,,] <- t(f1)
  
  # finding diagnosis label
  dx_group[i] <- phenotypic$DX_GROUP[which(phenotypic$SUB_ID == sub_id[i])]
}


# The first 73 are in group 1, the rest 98 are in group 2. 171 overall.
h.1 <- h[1:73,,]
h.2 <- h[74:171,,]

save(h.1, h.2, file="")
