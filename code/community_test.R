#libraries
library(igraph)
library(binom)

#global params
doALL<-FALSE
geneSet<- 'dry' #'cold' 

#data
dat<-read.table('../data/Feltus_etal_AdditionalFile3_EdgeLists.txt', header=TRUE)
if(geneSet=='dry'){
  genes<-read.csv('../data/DaveAllGxTDrought.csv')
}
if(geneSet=='cold'){
  genes<-read.csv('../data/HannahAllGxTCold.csv')
}

genes <- as.character(genes[,1])

net_genes <- c(as.character(dat[,2]),as.character(dat[,3]))
net_genes <- unique(net_genes)
not_in <- 1-length(which(! as.character(genes) %in% net_genes))/nrow(genes)

#making arab main_screen == repeat_screen graph
dat.graph<-dat[which(dat[,1]=='GLOBAL'),2:3]
dat.graph[,1:2]<-apply(dat.graph[,1:2],2,as.character)

gr<-graph.data.frame(dat.graph,directed=FALSE)

deg<-degree(gr)

cl.gr<-clusters(gr)

bet<-betweenness(gr)

close<-closeness(gr)

eigen<-evcent(gr)

dat.out<-data.frame(names(bet),deg,bet,close,eigen$vector)
dat.out[,1]<-as.character(dat.out[,1])
colnames(dat.out)<-c('gene','degree','betweenness','closeness','eigenvector')

coms <- leading.eigenvector.community(gr)
com_ids <- unique(coms$membership)

conf_ints <- binom.bayes(table(as.factor(coms$membership)),sum(table(as.factor(coms$membership))), conf.level=0.99) #we can do better than this

whole_prop <- log(as.numeric(table(as.factor(coms$membership))/sum(table(as.factor(coms$membership)))))
ge_prop <- as.numeric(table(as.factor(coms$membership)[which(coms$names %in% genes)])/sum(table(as.factor(coms$membership)[which(coms$names %in% genes)])))

com_out <- data.frame(read.csv("../data/gene_expression_net_feltus_main.csv"))
com_out_mt <- match(com_out$gene, coms$names)
com_out$community <- coms$membership[com_out_mt]

if(doALL == TRUE){
  com.res <- matrix(NA, ncol = 4, nrow=length(com_ids))
  colnames(com.res) <- c('degree','betweenness','closeness','eigenvector')
  
  counter <- 1
  for(c in com_ids){
    genes.c <- coms$names[which(coms$membership==c)]
    use.c <- which(dat.graph[,1] %in% genes.c & dat.graph[,2] %in% genes.c)
    gr.c <- graph.data.frame(dat.graph[use.c,],directed=FALSE)
    
    deg.c <- degree(gr.c)
    
    cl.gr.c <- clusters(gr.c)
    
    bet.c <- betweenness(gr.c)
    
    close.c <- closeness(gr.c)
    
    eigen.c <- evcent(gr.c)
    
    id.c <- rep(1,length(V(gr.c)$name))
    id.c[match(V(gr.c)$name, genes)] <- 2
    
    dat.test.c <- data.frame(names(bet.c),deg.c,bet.c,close.c,eigen.c$vector)
    dat.test.c[,1]<-as.character(dat.test.c[,1])
    colnames(dat.test.c)<-c('gene','degree','betweenness','closeness','eigenvector')
    
    #layout(matrix(1:4,ncol=2))
    for(i in 2:ncol(dat.test.c)){
      dat.i <- dat.test.c[which(is.finite(log(dat.test.c[,i]))==TRUE),]
      use <- which(dat.i[,1] %in% genes)
      id <- rep('nonGE',nrow(dat.i))
      id[use] <- 'GE'
      dat.i <- log(dat.i[,i])
      obs.i <- abs(mean(dat.i[use])-mean(dat.i[-use]))
      len1 <- c()
      for(j in 1:1000){
        use.j <- sample(1:length(dat.i),length(use),replace=TRUE)
        met.j <- abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
        len1 <- c(len1,met.j)
      }
      p.i <- round(length(which(len1>obs.i)),4)/j
      #boxplot(dat.i~id,main=paste0(colnames(dat.test.c)[i],' p value = ',p.i))
      if(p.i < 0.05){
        test.c.i <- by(dat.i, id, mean, na.rm=TRUE)
        diff.c.i <- test.c.i['GE'] - test.c.i['nonGE']
        com.res[counter, colnames(dat.test.c)[i]] <- diff.c.i
      }else{
        com.res[counter, colnames(dat.test.c)[i]] <- 0
      }	
    }	
    counter <- counter + 1
    
  }
  
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  barplot(com.res[which(com.res[,1]!=0),1], main = "degree", ylab = "GxE - Full Community")
  barplot(com.res[which(com.res[,2]!=0),2], main = "betweenness", ylab = "GxE - Full Community")
  barplot(com.res[which(com.res[,3]!=0),3], main = "closeness", ylab = "GxE - Full Community")
  barplot(com.res[which(com.res[,4]!=0),4], main = "eigenvector", ylab = "GxE - Full Community")
}