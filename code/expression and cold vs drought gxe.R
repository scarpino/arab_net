rm(list=ls())

#libraries
library(igraph)
library(glmnet)

#data
exp <- read.csv("../data/Expression Level/Yang&Gaut Table_S1.csv")

cold <- read.csv("../data/HannahAllGxTCold.csv")

cold_0.1 <- read.csv("../data/FDR0.1/Cold GxE FDR0.1.csv")


dry <- read.csv("../data/DaveAllGxTDrought.csv")

dry_0.1 <- read.csv("../data/FDR0.1/Drought GxE FDR0.1.csv")

IND <- rep(NA, nrow(exp))
IND[which(exp$ID %in% cold$Probe_Set_ID)] <- "cold"
IND[which(exp$ID %in% dry[,1])] <- "dry"

boxplot(log(exp$Expression.Level) ~ IND)

t.test(log(exp$Expression.Level) ~ IND)

metric <- diff(by(exp$Expression.Level, IND, mean))

sim_mets <- c()
for(i in 1:1000){
	IND.i <- rep(NA, nrow(exp))
	IND.i[sample(1:nrow(exp), length(which(IND == "cold")))] <- "cold"
	IND.i[sample(1:nrow(exp), length(which(IND == "dry")))] <- "dry"
	metric.i <- diff(by(exp$Expression.Level, IND.i, mean))
	sim_mets <- c(sim_mets, metric.i)
}

length(which(abs(sim_mets) > abs(metric)))/i

#re-sampling cold
prob <- exp$Expression.Level/by(exp$Expression.Level, IND, sum)[1]

use <- sample(which(IND == "cold"), size = 900, prob = 1/prob[which(IND == "cold")])

IND2 <- rep(NA, nrow(exp))
IND2[use] <- "cold"
IND2[which(exp$ID %in% dry[,1])] <- "dry"

t.test(log(exp$Expression.Level) ~ IND2)

dat.test<- read.csv("~/Dropbox (EmergentEpidemicsLab)/arab/arab shared/New Analyses 2.2017/gene_expression_net_feltus_main.csv")

genes <- as.character(exp$ID[which(IND2 == "cold")])

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