############
#Cold genes#
############

geneSet <- 'cold' 

#read in data
dat<-read.table('../data/Feltus_etal_AdditionalFile3_EdgeLists.txt', header=TRUE)
if(geneSet=='dry'){
	genes<-read.csv('../data/DaveAllGxTDrought.csv')
}
if(geneSet=='cold'){
	genes<-read.csv('../data/HannahAllGxTCold.csv')
}
genes<-as.character(genes[,1])


#making expression network
dat.graph <- dat[which(dat[,1] == 'GLOBAL'), 2:3]

gr <- graph.data.frame(dat.graph, directed=FALSE)

#calculating metrics
deg<-degree(gr)

cl.gr<-clusters(gr)

bet<-betweenness(gr)

close<-closeness(gr)

eigen<-evcent(gr)

id <- rep(1,length(V(gr)$name))
id[match(V(gr)$name, genes)] <- 2

cols <- c('gray','#b2182b')
size <- c(1,2)
V(gr)$color <- cols[id]
V(gr)$size <- size[id]

#stats
dat.test <-data.frame(names(bet),deg,bet,close,eigen$vector)
dat.test[,1]<-as.character(dat.test[,1])
colnames(dat.test)<-c('gene','degree','betweenness','closeness','eigenvector')

pvals <- list()
pvals[["cold"]] <- list()
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
	#boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
	#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
	#abline(v=obs.i,col='red',lty=3,lwd=3)
	pvals[["cold"]][[colnames(dat.test)[i]]] <- p.i
}

res.cold <- data.frame(names(bet), deg, bet, close, eigen$vector)
res.cold[,1] <- as.character(res.cold[,1])
colnames(res.cold) <- c('gene','degree','betweenness','closeness','eigenvector')
use<-which(res.cold[,"gene"] %in% genes)
id<-rep('nonGE',nrow(res.cold))
id[use]<-'GE'

res.cold <- data.frame(id, res.cold)

#80/20 testing
N <- 1000
draw.size <- floor(nrow(dat.test) * 0.8)
er.cold <- c()
gene.draw.cold <- list()
for(n in 1:N){
	gens.n <- sample(dat.test[,1], size = draw.size, replace = FALSE)
	in.n <- which(dat.test[,1] %in% gens.n)
	dat.n <- dat.test[in.n, ]
	GxE.n <- which(dat.n[,1] %in% genes)
	id.n <- rep(0, nrow(dat.n))
	id.n[GxE.n] <- 1
	mod.n <- glm(id.n~degree + eigenvector, family = binomial, data = dat.n)
	pred.n <- predict(mod.n, dat.test[-in.n,-1], type = "response")
	true.n <- rep(0, nrow(dat.test[-in.n,]))
	true.n[which(dat.test[-in.n,1] %in% genes)] <- 1
	er.n <- abs(pred.n-true.n)
	er.cold <- c(er.cold, length(which(er.n <= 0.05))/length(er.n))
	gene.draw.cold[[n]] <- dat.n[,1]
}

#community tests
coms <- leading.eigenvector.community(gr)
com_ids <- unique(coms$membership)

conf_ints <- binom.bayes(table(as.factor(coms$membership)),sum(table(as.factor(coms$membership))), conf.level=0.99) #we can do better than this

whole_prop <- as.numeric(table(as.factor(coms$membership))/sum(table(as.factor(coms$membership))))
ge_prop <- as.numeric(table(as.factor(coms$membership)[which(coms$names %in% genes)])/sum(table(as.factor(coms$membership)[which(coms$names %in% genes)])))
ge_cols <- rep("#878787", length(ge_prop))
sig_diff <- which(ge_prop > conf_ints[,'upper']|ge_prop < conf_ints[,'lower'])
#ge_cols[sig_diff] <- "#b2182b"
sigs <- rep(0, length(ge_prop))
sigs[sig_diff] <- 1

com.res <- matrix(NA, ncol = 4, nrow=length(com_ids))
colnames(com.res) <- c('degree','betweenness','closeness','eigenvector')

dat.props <- data.frame(whole_prop, as.numeric(table(as.factor(coms$membership))), ge_prop, as.numeric(table(as.factor(coms$membership)[which(coms$names %in% genes)])), sigs)
colnames(dat.props) <- c("full_prop", "full_raw", "cold_prop", "cold_raw", "sig_cold")

diffs <- dat.props$full_prop-dat.props$cold_prop
use.coms <- which(abs(diffs) > cut_sig & dat.props$sig_cold == 1)

for(c in use.coms){
	dat.use <- which(coms$membership == c)
	dat.out <- data.frame(coms$names[dat.use])
	if(diffs[c] > 0){
		out.name <- "overrepresented"
	}else{
		out.name <- "underrepresented"
	}
	write.csv(dat.out, file = paste0(c,"_cold_",out.name,"_", cut_sig, ".csv"), row.names = FALSE, quote = FALSE)
}

if(do_com_test == TRUE){
	counter <- 1
	big_coms_res_cold <- list()
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
			if(p.i < 1){
				test.c.i <- by(exp(dat.i), id, median, na.rm=TRUE)
				diff.c.i <- test.c.i['GE'] - test.c.i['nonGE']
				com.res[counter, colnames(dat.test.c)[i]] <- diff.c.i
			}else{
				com.res[counter, colnames(dat.test.c)[i]] <- 0
			}	
		}	
		counter <- counter + 1
	
	}
	big_coms_res_cold[[c]] <- com.res
	com.est.cold <- list()
	for(col in 1:ncol(com.res)){
		com.est.cold[[colnames(com.res)[col]]] <- sum(com.res[,col] * ge_prop, na.rm = TRUE)
	}
}

#for graph compare
genes.cold <- coms$names[which(coms$membership==which.max(ge_prop))]
use.cold <- which(dat.graph[,1] %in% genes.cold & dat.graph[,2] %in% genes.cold)

gr.cold.com <- graph.data.frame(dat.graph[use.cold,], directed = FALSE)

sizes <- c(1.5,3)
cols.com.plot <- c("#bababa","#b2182b","#2166ac") #plain, dry, cold

colors.cold <- rep(1, length(get.vertex.attribute(gr.cold.com, "name")))
sizes.cold <- rep(1, length(get.vertex.attribute(gr.cold.com, "name")))

colors.cold[ which(get.vertex.attribute(gr.cold.com, "name") %in% genes)] <- 3
sizes.cold[ which(get.vertex.attribute(gr.cold.com, "name") %in% genes)] <- 2


###############
#Drought genes#
###############

geneSet <- 'dry' 

#read in data
if(geneSet=='dry'){
	genes<-read.csv('../data/DaveAllGxTDrought.csv')
}
if(geneSet=='cold'){
	genes<-read.csv('../data/HannahAllGxTCold.csv')
}
genes<-as.character(genes[,1])


#making expression network
dat.graph <- dat[which(dat[,1] == 'GLOBAL'), 2:3]

gr <- graph.data.frame(dat.graph, directed=FALSE)

#calculating metrics
deg<-degree(gr)

cl.gr<-clusters(gr)

bet<-betweenness(gr)

close<-closeness(gr)

eigen<-evcent(gr)

id <- rep(1,length(V(gr)$name))
id[match(V(gr)$name, genes)] <- 2

cols <- c('gray','#b2182b')
size <- c(1,2)
V(gr)$color <- cols[id]
V(gr)$size <- size[id]

#stats
dat.test <-data.frame(names(bet),deg,bet,close,eigen$vector)
dat.test[,1]<-as.character(dat.test[,1])
colnames(dat.test)<-c('gene','degree','betweenness','closeness','eigenvector')

pvals[["dry"]] <- list()
for(i in 2:ncol(dat.test)){
	dat.i<-dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
	use<-which(dat.i[,1] %in% genes)
	id<-rep('nonGE',nrow(dat.i))
	id[use]<-'GE'
	dat.i<-dat.i[,i]
	obs.i<-abs(mean(dat.i[use])-mean(dat.i[-use]))
	len1<-c()
	for(j in 1:10000){
		use.j<-sample(1:length(dat.i),length(use),replace=TRUE)
		met.j<-abs(mean(dat.i[use.j])-mean(dat.i[-use.j]))
		len1<-c(len1,met.j)
	}
	p.i<-round(length(which(len1>obs.i)),4)/j
	#boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
	#hist(len1,ain=paste0(colnames(dat.test)[i],' p value = ',p.i),col='gray')
	#abline(v=obs.i,col='red',lty=3,lwd=3)
	pvals[["dry"]][[colnames(dat.test)[i]]] <- p.i
}

res.dry <- data.frame(names(bet), deg, bet, close, eigen$vector)
res.dry[,1] <- as.character(res.dry[,1])
colnames(res.dry) <- c('gene','degree','betweenness','closeness','eigenvector')
use<-which(res.dry[,"gene"] %in% genes)
id<-rep('nonGE',nrow(res.dry))
id[use]<-'GE'

res.dry <- data.frame(id, res.dry)

#80/20 testing
N <- 1000
draw.size <- floor(nrow(dat.test) * 0.8)
er.dry <- c()
gene.draw.dry <- list()
for(n in 1:N){
	gens.n <- sample(dat.test[,1], size = draw.size, replace = FALSE)
	in.n <- which(dat.test[,1] %in% gens.n)
	dat.n <- dat.test[in.n, ]
	GxE.n <- which(dat.n[,1] %in% genes)
	id.n <- rep(0, nrow(dat.n))
	id.n[GxE.n] <- 1
	mod.n <- glm(id.n~degree + eigenvector, family = binomial, data = dat.n)
	pred.n <- predict(mod.n, dat.test[-in.n,-1], type = "response")
	true.n <- rep(0, nrow(dat.test[-in.n,]))
	true.n[which(dat.test[-in.n,1] %in% genes)] <- 1
	er.n <- abs(pred.n-true.n)
	er.dry <- c(er.dry, length(which(er.n <= 0.05))/length(er.n))
	gene.draw.dry[[n]] <- dat.n[,1]
}


#community tests - dry
ge_prop <- as.numeric(table(as.factor(coms$membership)[which(coms$names %in% genes)])/sum(table(as.factor(coms$membership)[which(coms$names %in% genes)])))
ge_cols <- rep("#878787", length(ge_prop))
sig_diff <- which(ge_prop > conf_ints[,'upper']|ge_prop < conf_ints[,'lower'])
#ge_cols[sig_diff] <- "#b2182b"

sigs <- rep(0, length(ge_prop))
sigs[sig_diff] <- 1

dat.props <- data.frame(dat.props, ge_prop, as.numeric(table(as.factor(coms$membership)[which(coms$names %in% genes)])), sigs)
colnames(dat.props) <- c("full_prop", "full_raw", "cold_prop", "cold_raw", "sig_cold", "dry_prop", "dry_raw", "dry_sig")
#dat.props[,"full_prop"] <- exp(dat.props[,"full_prop"])

com.res <- matrix(NA, ncol = 4, nrow=length(com_ids))
colnames(com.res) <- c('degree','betweenness','closeness','eigenvector')

diffs <- dat.props$full_prop-dat.props$dry_prop
use.coms <- which(abs(diffs) > cut_sig & dat.props$dry_sig == 1)

for(c in use.coms){
	dat.use <- which(coms$membership == c)
	dat.out <- data.frame(coms$names[dat.use])
	if(diffs[c] > 0){
		out.name <- "overrepresented"
	}else{
		out.name <- "underrepresented"
	}
	write.csv(dat.out, file = paste0(c,"_dry_",out.name,"_", cut_sig,".csv"), row.names = FALSE, quote = FALSE)
}

if(do_com_test == TRUE){
	big_coms_res_dry <- list()
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
			if(p.i < 1){
				test.c.i <- by(exp(dat.i), id, median, na.rm=TRUE)
				diff.c.i <- test.c.i['GE'] - test.c.i['nonGE']
				com.res[counter, colnames(dat.test.c)[i]] <- diff.c.i
			}else{
				com.res[counter, colnames(dat.test.c)[i]] <- 0
			}	
		}	
		counter <- counter + 1
	}
	big_coms_res_dry[[c]] <- com.res
	com.est.dry <- list()
	for(col in 1:ncol(com.res)){
		com.est.dry[[colnames(com.res)[col]]] <- sum(com.res[,col] * ge_prop, na.rm = TRUE)
	}
}

#for graph compare
genes.dry <- coms$names[which(coms$membership==which.max(ge_prop))]
use.dry <- which(dat.graph[,1] %in% genes.dry & dat.graph[,2] %in% genes.dry)

gr.dry.com <- graph.data.frame(dat.graph[use.dry,], directed = FALSE)

colors.dry <- rep(1, length(get.vertex.attribute(gr.dry.com, "name")))
sizes.dry <- rep(1, length(get.vertex.attribute(gr.dry.com, "name")))

colors.dry[ which(get.vertex.attribute(gr.dry.com, "name") %in% genes)] <- 2
sizes.dry[ which(get.vertex.attribute(gr.dry.com, "name") %in% genes)] <- 2