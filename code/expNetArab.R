rm(list=ls())

#libraries
library(igraph)

#global params
doALL<-FALSE
geneSet<- 'cold' #'dry' 

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
write.csv(dat.out,file='geneExpressionNet.csv',row.names=FALSE)

null_set <- read.csv('10_2011_freezing_ath1_AllAcc.csv')[,1]
null_set <-  unlist(strsplit(as.character(null_set), '_'))
null_set <- null_set[seq(1,length(null_set), by = 2)]

#stats
#dat.test<-dat.out[-which(dat.out$betweenness==0),]
dat.test<-dat.out
layout(matrix(1:4,ncol=2))
for(i in 2:ncol(dat.test)){
	dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
	use<-which(dat.i[,1] %in% genes)
	id<-rep('nonGE',nrow(dat.i))
	id[use]<-'GE'
	null_set.i <- which(dat.i[,1] %in% null_set)
	dat.i<-log(dat.i[,i])
	obs.i<-abs(mean(dat.i[use])-mean(dat.i[-use]))
	len1<-c()
	for(j in 1:10000){
		use.j<-sample(null_set.i,length(use),replace=TRUE)
		met.j<-abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
		len1<-c(len1,met.j)
	}
	p.i<-round(length(which(len1>obs.i)),4)/j
	boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
	#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
	#abline(v=obs.i,col='red',lty=3,lwd=3)
}


if(doALL==TRUE){
	nets<-unique(dat[,1])
	setwd(paste0('allNets/', geneSet,'/'))
	P<-matrix(NA,ncol=4,nrow=length(nets))
	E<-matrix(NA,ncol=4,nrow=length(nets))
	colnames(P)<-c('degree','betweenness','closeness','eigenvector')
	colnames(E)<-c('degree','betweenness','closeness','eigenvector')
	V_ge <- P
	V_non_ge <- P
	n_ge <- P
	for(n in 1:length(nets)){
		dat.graph<-dat[which(dat[,1]==nets[n]),2:3]
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
	write.csv(dat.out,file=paste0(nets[n],' geneExpressionNet.csv'),row.names=FALSE)
		
		#stats
		#dat.test<-dat.out[-which(dat.out$betweenness==0),]
		#pdf(paste0(nets[n],'.pdf'))
		#layout(matrix(1:4,ncol=2))
	  dat.test<-dat.out
		for(i in 2:ncol(dat.test)){
			dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
			use<-which(dat.i[,1] %in% genes)
			id<-rep('nonGE',nrow(dat.i))
			id[use]<-'GE'
			dat.i<-log(dat.i[,i])
			obs.i<-abs(mean(dat.i[use])-mean(dat.i[-use]))
			len1<-c()
			do_iter <- TRUE
			if(do_iter == TRUE){
			  for(j in 1:1000){
			    use.j<-sample(1:length(dat.i),length(use),replace=TRUE)
			    met.j<-abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
			    len1<-c(len1,met.j)
			  }
			  p.i<-round(length(which(len1>obs.i)),4)/j
			  boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
			  P[n,(i-1)]<-p.i
			}
			
			E[n,(i-1)]<-mean(dat.i[use])-mean(dat.i[-use])
			V_ge[n,(i-1)] <- mean(dat.i[use])
			V_non_ge[n,(i-1)] <- mean(dat.i[-use])
			n_ge[n,(i-1)] <- length(use)
			#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
			#abline(v=obs.i,col='red',lty=3,lwd=3)
		}#end for i
		#dev.off()
	}#end for n
	dat.out <- data.frame(as.character(nets),E,P, V_ge, V_non_ge, n_ge)
	write.csv(dat.out,file=paste0(geneSet,'_GILS.csv'))
	
	w <- n_ge[-1,"eigenvector"]/sum(n_ge[-1,"eigenvector"])
	sum(V_ge[-1,"eigenvector"]*w)
	w_gil <- table(dat$NETWORK)[-length(table(dat$NETWORK))]/sum(table(dat$NETWORK)[-length(table(dat$NETWORK))])
	sum(V_non_ge[-1,"eigenvector"]*w_gil)
	
	w <- n_ge[-1,"degree"]/sum(n_ge[-1,"degree"])
	sum(V_ge[-1,"degree"]*w)
	w_gil <- table(dat$NETWORK)[-length(table(dat$NETWORK))]/sum(table(dat$NETWORK)[-length(table(dat$NETWORK))])
	sum(V_non_ge[-1,"degree"]*w_gil)
	
	gils<-c(4, 6, 12, 18, 19, 29,40,46,47,51,55, 59,66,71,75,77,80,82,86) #stress gils
	#gils.leaf = c(6,12,18,19,29, 40, 46,47,51,55,59,75,77,80,82,86)
	gils<-paste0('GIL',gils)
	use.e<-which(nets%in%gils)
	pdf('gilSummary.pdf')
	layout(matrix(1:4,ncol=2))
	for(e in 1:ncol(E)){
		hist(E[use.e,e],breaks=10,main=paste0(colnames(P)[e],' prop sig = ',round(length(E[use.e,e][which(P[use.e,e]<0.05)])/length(E[use.e,e]),3)),col='#00000050',xlab='mean change GE - nonGE')
		hist(E[use.e,e][which(P[use.e,e]<0.05)],add=TRUE,col='#000000',breaks=10)
	}#end for e
	dev.off()
}#end if doAll