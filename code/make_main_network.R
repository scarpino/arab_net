#libraries
library(igraph)

#data
dat <- read.table('../data/Feltus_etal_AdditionalFile3_EdgeLists.txt', header=TRUE)

#making arab graph
dat.graph <- dat[which(dat[,1]=='GLOBAL'),2:3]
dat.graph[,1:2] <- apply(dat.graph[,1:2], 2, as.character)

gr <- graph.data.frame(dat.graph,directed=FALSE)

deg <- degree(gr)

cl.gr <- clusters(gr)

bet <- betweenness(gr)

close <- closeness(gr)

eigen <- evcent(gr)

dat.out <- data.frame(names(bet),deg,bet,close,eigen$vector)
dat.out[,1] <- as.character(dat.out[,1])
colnames(dat.out) <- c('gene','degree','betweenness','closeness','eigenvector')
#write.csv(dat.out,file='gene_expression_net_feltus_main.csv',row.names=FALSE)