#libraries
library(igraph)

#global params
doALL<-FALSE
geneSet<- 'cold' #'dry' # #

#data
dat<-read.table('../data/Feltus_etal_AdditionalFile3_EdgeLists.txt', header=TRUE)
if(geneSet=='dry'){
	genes<-read.csv('../data/DaveAllTDrought.csv')
}
if(geneSet=='cold'){
	genes<-read.csv('../data/HannahAllTCold.csv')
}
genes<-as.character(genes[,1])

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

#stats
dat.test<-dat.out
layout(matrix(1:4,ncol=2))
for(i in 2:ncol(dat.test)){
	dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
	use<-which(dat.i[,1] %in% genes)
	id<-rep('nonE',nrow(dat.i))
	id[use]<-'E'
	null_set.i <- which(dat.i[,1] %in% genes)
	dat.i<-log(dat.i[,i])
	#dat.i<-dat.i[,i]
	obs.i<-abs(mean(dat.i[use])-mean(dat.i[-use]))
	len1<-c()
	for(j in 1:1000){
		use.j<-sample(use,length(use),replace=TRUE)
		met.j<-abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
		len1<-c(len1,met.j)
	}
	p.i<-round(length(which(len1>obs.i)),4)/j
	boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
	#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
	#abline(v=obs.i,col='red',lty=3,lwd=3)
}

#overlap of dry and cold genes
gd <- as.character(read.csv('../data/DaveAllTDrought.csv')[,1])
gc <- as.character(read.csv('../data/HannahAllTCold.csv')[,1])

use <- which(gd %in% gc)

g_overlap <- gd[use]

#stats
use_gd <- which(dat.out[,1] %in% gd)
use_gc <- which(dat.out[,1] %in% gc)
use_tot <- unique(c(use_gd, use_gc))

dat.test <- dat.out[use_tot,]
layout(matrix(1:4,ncol=2))
for(i in 2:ncol(dat.test)){
	dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
	use<-which(dat.i[,1] %in% g_overlap)
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
	boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i), range = 0)
	#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
	#abline(v=obs.i,col='red',lty=3,lwd=3)
}




#E vs. GxE
gxe_cold <- as.character(read.csv('../data/HannahAllGxTCold.csv')[,1])
use_gd <- which(dat.out[,1] %in% gd)
use_gc <- which(dat.out[,1] %in% gc)
use_gxe_c <- which(dat.out[,1] %in% gxe_cold)
use_tot <- unique(c(use_gxe_c, use_gc,  use_gd))

coms <- fastgreedy.community(graph = gr)
use_cold_ge <- coms$membership[which(coms$names %in% gxe_cold)]

layout(matrix(1:2, nrow = 1))
barplot(table(use_cold_ge)[order(table(use_cold_ge), decreasing = TRUE)]/sum(table(use_cold_ge)), ylim = c(0,0.5), main = "Cold eGxE", xlab = "Community ID", ylab = "Prop.")

use_cold_e <- coms$membership[which(coms$names %in% gc)]
barplot(table(use_cold_e)[order(table(use_cold_e), decreasing = TRUE)]/sum(table(use_cold_e)), ylim = c(0, 0.5), main = "Cold eE", xlab = "Community ID", ylab = "Prop.")

p_e <- table(use_cold_e)/sum(table(use_cold_e))
p_ge <- table(use_cold_ge)/sum(table(use_cold_ge))
ent.diff <- entropy(p_e) - entropy(p_ge)

ents <- c()
for(i in 1:1000){
  e.i <- sample(coms$membership, length(use_cold_e))
  ge.i <- sample(coms$membership, length(use_cold_ge))
  
  p_e.i <- table(e.i)/sum(table(e.i))
  p_ge.i <- table(ge.i)/sum(table(ge.i))
  
  ent.diff.i <- entropy(p_e.i) - entropy(p_ge.i)
  
  ents <- c(ents, ent.diff.i)
}

length(which(abs(ents) > abs(ent.diff)))/1000

dat.test <- dat.out[use_tot,]
layout(matrix(1:4,ncol=2))
for(i in 2:ncol(dat.test)){
	dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
	use<-which(dat.i[,1] %in% gxe_cold)
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
	boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i), range = 0)
	#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
	#abline(v=obs.i,col='red',lty=3,lwd=3)
}


	
nets<-unique(dat[,1])
if(doALL==TRUE){
	setwd(paste0('allNets/', geneSet,'/'))
	P<-matrix(NA,ncol=4,nrow=length(nets))
	E<-matrix(NA,ncol=4,nrow=length(nets))
	colnames(P)<-c('degree','betweenness','closeness','eigenvector')
	colnames(E)<-c('degree','betweenness','closeness','eigenvector')
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
	#write.csv(dat.out,file=paste0(nets[n],' geneExpressionNet.csv'),row.names=FALSE)
		
		#stats
		dat.test<-dat.out[-which(dat.out$betweenness==0),]
		pdf(paste0(nets[n],'.pdf'))
		layout(matrix(1:4,ncol=2))
		for(i in 2:ncol(dat.test)){
			dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
			use<-which(dat.i[,1] %in% genes)
			id<-rep('nonGE',nrow(dat.i))
			id[use]<-'GE'
			dat.i<-log(dat.i[,i])
			obs.i<-abs(mean(dat.i[use])-mean(dat.i[-use]))
			len1<-c()
			for(j in 1:10000){
				use.j<-sample(1:length(dat.i),length(use),replace=TRUE)
				met.j<-abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
				len1<-c(len1,met.j)
			}
			p.i<-round(length(which(len1>obs.i)),4)/j
			boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
			P[n,(i-1)]<-p.i
			E[n,(i-1)]<-mean(dat.i[use])-mean(dat.i[-use])
			#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
			#abline(v=obs.i,col='red',lty=3,lwd=3)
		}#end for i
		dev.off()
	}#end for n
	#write.csv(P,file=paste0(nets[n],'_P.pdf'))
	#write.csv(E,file=paste0(nets[n],'_E.pdf'))
	
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