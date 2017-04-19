#libraries
library(igraph)
library(binom)

#running analyses
sizes <- c(1.5,3)
do_com_test <- FALSE
cut_sig <- 0.03 #difference in gene proportions for sig. communities.

source("arab acc funcs.R")
source("arab_overunder_funcs.R")
#Figure 1
quartz(width = 5, height = 5)
par(mar=c(2,4,1,0.2))
layout(matrix(1:4, ncol = 2, byrow = TRUE))
boxplot(log(res.dry[,"degree"])~res.dry[,"id"], col = c("#b2182b","#f4a582"), yaxt = "n", ylab = "Degree", names = c(" ", " "), ylim = c(0,7), range = 0)
axis(2,at=seq(0, 7, length.out = 6),labels=round(exp(seq(0, 7, length.out = 6)),0),las = 2)
#text(1.5,6.4,"Dry")
title(main = "Dry")

par(mar=c(2,2,1,2))
boxplot(log(res.cold[,"degree"])~res.cold[,"id"], col = c("#2166ac","#92c5de"), yaxt = "n", names = c(" ", " "), ylim = c(0,7), range = 0)
axis(2,at=seq(0, 7, length.out = 6),labels=rep("", 6),las = 2)
#text(1.5,6.4,"Cold")
title(main = "Cold")

par(mar=c(2,4,1,0.2))
use.eigen <- which(is.finite(log(res.dry[,"eigenvector"]))==TRUE)
boxplot(log(res.dry[use.eigen,"eigenvector"])~res.dry[use.eigen,"id"], col = c("#b2182b","#f4a582"), yaxt = "n", ylab = "Eigenvector centrality (log)", names = c("G x E", "non G x E"), ylim = c(-43, 0), range = 0)
axis(2,at=seq(-43, 0, length.out = 6),labels=round(seq(-43, 0, length.out = 6),0),las = 2)

par(mar=c(2,2,1,2))
use.eigen <- which(is.finite(log(res.cold[,"eigenvector"]))==TRUE)
boxplot(log(res.cold[use.eigen,"eigenvector"])~res.cold[use.eigen,"id"], col = c("#2166ac","#92c5de"), yaxt = "n", names = c("G x E", "non G x E"), pch = 16, ylim = c(-43, 0), range = 0)
axis(2,at=seq(-43, 0, length.out = 6),labels=rep("", 6),las = 2)


#Figure 2
png("dry_comm.png", width = 1800, height = 1800)
par(mar = c(0,1,4,1))
plot(gr.dry.com, vertex.label = NA, vertex.size = sizes[sizes.dry], vertex.color = cols.com.plot[colors.dry], edge.color = "#e0e0e095", edge.width = 1, main = "Principle Dry Community", vertex.frame.color = NA, layout = layout_nicely(gr.dry.com))
dev.off()


png("cold_comm.png", width = 1800, height = 1800)
par(mar = c(0,1,4,1))
plot(gr.cold.com,  vertex.label = NA, vertex.size = sizes[sizes.cold], vertex.color = cols.com.plot[colors.cold], edge.color = "#e0e0e095", edge.width = 1, main = "Principle Cold Community", vertex.frame.color = NA, layout = layout_nicely(gr.cold.com))
dev.off()

pdf("F2.pdf")
#quartz(width = 10, height = 4)
layout(matrix(1:2, ncol = 2))
par(mar = c(0,0.1,4,0.5))
plot(gr.dry.com, vertex.label = NA, vertex.size = sizes[sizes.dry], vertex.color = cols.com.plot[colors.dry], edge.color = "#e0e0e055", edge.width = 0.5, main = "Principle Dry Community", vertex.frame.color = NA, layout = layout_nicely(gr.dry.com))

par(mar = c(0,0.5,4,0.1))
plot(gr.cold.com,  vertex.label = NA, vertex.size = sizes[sizes.cold], vertex.color = cols.com.plot[colors.cold], edge.color = "#e0e0e055", edge.width = 0.5, main = "Principle Cold Community", vertex.frame.color = NA, layout = layout_nicely(gr.cold.com))
dev.off()

#testing with TPC2015-00910
clust_TPC <- read.csv("../data/TPC2015_00910/clusters_TPC2015_00910.csv")[,1:4]
length(which(genes.dry %in% clust_TPC$Locus.ID))

ALG22 <- "AT2G22540"
gen.unq <- unique(c(as.character(dat[,2]), as.character(dat[,3])))
scores.dry <- rep(NA, length(gen.unq))
for(i in 1:length(scores.dry)){
	use.i <- which(dat[,2] == gen.unq[i] | dat[,3] == gen.unq[i])
	rm.i <- which(genes.dry == gen.unq[i])
	if(length(rm.i) > 0){
		genes.dry.test <- genes.dry[-rm.i]
	}else{
		genes.dry.test <- genes.dry
	}
	score.i <- length(which(unique(c(as.character(dat[use.i,2]), as.character(dat[use.i,3]))) %in% genes.dry.test))
	scores.dry[i] <- score.i
}

scores.cold <- rep(NA, length(gen.unq))
for(i in 1:length(scores.cold)){
	use.i <- which(dat[,2] == gen.unq[i] | dat[,3] == gen.unq[i])
	rm.i <- which(genes.cold == gen.unq[i])
	if(length(rm.i) > 0){
		genes.cold.test <- genes.cold[-rm.i]
	}else{
		genes.cold.test <- genes.cold
	}
	score.i <- length(which(unique(c(as.character(dat[use.i,2]), as.character(dat[use.i,3]))) %in% genes.cold.test))
	scores.cold[i] <- score.i
}

ord.dry <- order(scores.dry, decreasing = TRUE)
ord.cold <- order(scores.cold, decreasing = TRUE)
out.scores <- data.frame(gen.unq[ord.dry[1:10]], scores.dry[ord.dry[1:10]], gen.unq[ord.cold[1:10]], scores.cold[ord.cold[1:10]])
colnames(out.scores) <- c("Dry.hits", "Dry.Connex", "Cold.hits", "Cold.connex")
#write.csv(out.scores, file = "gene_connex_hits.csv", row.names = FALSE, quote = FALSE)

layout(matrix(1:2, nrow = 1))
hist(log(scores.dry), xaxt = "n", xlab = "Number of connects to drought GxE", freq = F, ylim = c(0,1), main = "Drought GxE", col = cols.com.plot[2])
at.x <- seq(0, max(log(scores.dry), na.rm = TRUE), length.out = 5)
axis(1, at = at.x, labels = round(exp(at.x), 0))

hist(log(scores.cold), xaxt = "n", xlab = "Number of connects to cold GxE", freq = F, ylim = c(0,1), main = "Cold GxE", col = cols.com.plot[3])
at.x <- seq(0, max(log(scores.cold), na.rm = TRUE), length.out = 5)
axis(1, at = at.x, labels = round(exp(at.x), 0))

#Figure 3
overunder <- ggplot(commsFilter)+
  geom_hline(yintercept=0, linetype = 2)+
  geom_point(aes(y=relDiff,x=size, color=colCat), size=2)+
  scale_x_log10(limits=c(10,1000))+
  facet_wrap(~treatment)+theme_few(base_size = 16)+
  theme(legend.position="none")+
  scale_color_manual(values=c( "darkorange2", "grey90","chartreuse4"))+
  labs(x="Community size", y="Relative representation of eGxE genes")

cowplot::ggdraw()+
  cowplot::draw_plot(overunder)+
  cowplot::draw_text("over", x=0.45, y=0.85, color="chartreuse4", size = 11)+
  cowplot::draw_text("under", x=0.45, y=0.4, color="darkorange2", size = 11)
ggsave("F3.pdf", height = 5, width=7)
