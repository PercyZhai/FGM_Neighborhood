##################
#  ADHD Data Reader
#  Reads in NYU data and generate n*p*T array h[i,j,t]
##################
library(stringr)

file.path <- "Desktop/GroupLasso/ADHD200_AAL_TCs_filtfix/NYU"
#file.list <- list.files(file.path, "NYU?")
phenotypic <- read.csv(paste(file.path, "NYU_phenotypic.csv", sep="/"))

# extract ScanDir ID
scandir.id <- phenotypic[which(phenotypic$QC_Rest_1 == "1") ,1]
dx.group <- phenotypic[which(phenotypic$QC_Rest_1 == "1") ,6]
scandir.id <- str_pad(scandir.id, 7, pad="0")

n <- length(scandir.id)
p <- 116
tau <- 172

h <- array(0, c(n,p,tau))
dx_group <- rep(0,n)

for(i in 1:n){
  id <- scandir.id[i]
  f1 <- read.delim(paste(file.path, "/", id, "/snwmrda", 
                         id,"_session_1_rest_1_aal_TCs.1D", sep=""))
  f1 <- as.matrix(f1[ , -c(1,2)])
  # print(dim(f1)) everyone 172*116
  h[i,,] <- t(f1)
}

# ADHD has DX!=0, control has DX==0
h.1 <- h[which(dx.group != "0"),,] # ADHD
h.2 <- h[which(dx.group == "0"),,] # ctrl

print(dim(h.1))
print(dim(h.2))

save(h.1, h.2, file="Desktop/GroupLasso/Real data analysis/Data/time.series.ADHD.RData")
