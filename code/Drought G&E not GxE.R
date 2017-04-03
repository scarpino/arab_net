#libraries
library(igraph)

#data
dat.test<- read.csv("../data/gene_expression_net_feltus_main.csv")

genes <- read.csv("../data/G&E/Drought G&E NOT GxE.csv")

genes <- as.character(genes[,1])

layout(matrix(1:4,ncol=2))
for(i in 2:ncol(dat.test)){
  dat.i <- dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
  use <- which(dat.i[,1] %in% genes)
  id <- rep('nonGE',nrow(dat.i))
  id[use] <- 'GE'
  dat.i <- log(dat.i[,i])
  obs.i <- abs(mean(dat.i[use])-mean(dat.i[-use]))
  len1 <- c()
  for(j in 1:1000){
    use.j <- sample(1:length(id), length(use), replace=TRUE)
    met.j <- abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
    len1 <- c(len1,met.j)
  }
  p.i <- round(length(which(len1>obs.i)),4)/j
  boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
}
