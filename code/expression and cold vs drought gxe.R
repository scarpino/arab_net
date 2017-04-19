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
dat.test<- read.csv("../data/gene_expression_net_feltus_main.csv")

t_test_stats <- rep(NA, 1000)
pval_deg <- rep(NA, 1000)
pval_eigen <- rep(NA, 1000)
diff_deg <- rep(NA, 1000)
diff_eigen <- rep(NA, 1000)
for(n in 1:1000){

  use <- sample(which(IND == "cold"), size = 750, prob = 1/prob[which(IND == "cold")])
  
  IND2 <- rep(NA, nrow(exp))
  IND2[use] <- "cold"
  IND2[which(exp$ID %in% dry[,1])] <- "dry"
  
  test.i <- t.test(log(exp$Expression.Level) ~ IND2)
  t_test_stats[n] <- test.i$p.value
  
  genes <- as.character(exp$ID[which(IND2 == "cold")])
  
  layout(matrix(1:4,ncol=2))
  for(i in c(2,5)){
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
    if(i == 2){
      pval_deg[n] <- p.i
      diff_deg[n] <- met.j
    }
    if(i == 5){
      pval_eigen[n] <- p.i
      diff_eigen[n] <- met.j
    }
    #boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
  }
}

summary(t_test_stats)
summary(exp(diff_deg))
summary(diff_eigen)
summary(pval_eigen)
summary(pval_deg)
1 - length(which(pval_eigen < 0.05))/n
1 - length(which(pval_deg < 0.05))/n