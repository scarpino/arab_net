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

#testing for expression
dat.test<- read.csv("../data/gene_expression_net_feltus_main.csv")
genes <- as.character(exp$ID[which(IND == "cold")])

#logistic regression
rep_0 <- function(x){
  x[which(x == 0)] <- NA
  return(x)
}

dat.reg <- dat.test
dat.reg$expression <- exp$Expression.Level[match(dat.reg$gene, exp$ID)]
dat.reg <- dat.reg[,-1]

dat.reg <- apply(dat.reg, 2, rep_0)
dat.reg <- log(dat.reg)
GxE_ind <- rep(0, nrow(dat.reg)) 
GxE_ind[which(dat.test$gene %in% genes)] <- 1

m <- glm(GxE_ind~-1 + dat.reg)
summary(m)

#write.csv(summary(m)$coefficients, file = "expression_regression.csv")

#resampling
run_resample <- FALSE
if(run_resample == TRUE){
  t_test_stats <- rep(NA, 1000)
  pval_deg <- rep(NA, 1000)
  pval_eigen <- rep(NA, 1000)
  ge_deg <- rep(NA, 1000)
  ge_eigen <- rep(NA, 1000)
  non_ge_deg <- rep(NA, 1000)
  non_ge_eigen <- rep(NA, 1000)
  
  cold_exp <- exp$Expression.Level[which(IND == "cold")]
  cold_prob <- 1/log(cold_exp)
  cold_prob <- cold_prob/sum(cold_prob)
  pb <- txtProgressBar(1, 1000, style=3)
  for(n in 1:1000){
    
    use <- sample(which(IND == "cold"), prob = cold_prob, size = 0.5*length(which(IND =="cold")), replace = FALSE)
    
    IND2 <- rep(NA, nrow(exp))
    IND2[use] <- "cold"
    IND2[which(exp$ID %in% dry[,1])] <- "dry"
    genes_n <- as.character(exp$ID[which(IND2 == "cold")])
    
    test.i <- t.test(log(exp$Expression.Level) ~ IND2)
    t_test_stats[n] <- test.i$statistic
    
    layout(matrix(1:4,ncol=2))
    for(i in c(2,5)){
      dat.i <- dat.test[which(is.finite(log(dat.test[,i]))==TRUE),]
      rm <- which(!genes %in% genes_n)
      rm_out <- which(dat.i[,1] %in% rm)
      dat.i <- dat.i[-rm, ]
      use <- which(dat.i[,1] %in% genes_n)
      id <- rep('nonGE',nrow(dat.i))
      id[use] <- 'GE'
      dat.i <- log(dat.i[,i])
      #dat.i <- dat.i[,i]
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
        non_ge_deg[n] <- exp(mean(dat.i[-use.j]))
        ge_deg[n] <- exp(mean(dat.i[use.j]))
      }
      if(i == 5){
        pval_eigen[n] <- p.i
        non_ge_eigen[n] <- exp(mean(dat.i[-use.j]))
        ge_eigen[n] <- exp(mean(dat.i[use.j]))
      }
      #boxplot(dat.i~id,main=paste0(colnames(dat.test)[i],' p value = ',p.i))
    }
    setTxtProgressBar(pb, n)
  }
  
  quartz(width = 5, height = 5)
  par(mar=c(2,4,1,0.2))
  layout(matrix(1:2, ncol = 1))
  boxplot(log(ge_deg),log(non_ge_deg), col = c("#2166ac","#92c5de"), ylab = "Degree (log)", names = c("eGxE", "non eGxE"), range = 0)
  title(main = "Cold")
  boxplot(log(ge_eigen),log(non_ge_eigen), col = c("#2166ac","#92c5de"), ylab = "Eigenvector centrality (log)", names = c("eGxE", "non eGxE"), range = 0)
}