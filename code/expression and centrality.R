rm(list=ls())

#libraries
library(igraph)
library(glmnet)

#data
dat.test<- read.csv("../data/gene_expression_net_feltus_main.csv")

exp <- read.csv("../data/Expression Level/Yang&Gaut Table_S1.csv")

mt <- match(exp$ID, dat.test$gene)

x1 <- exp$Ka.Ks
x2 <- exp$Expression.Level
x3 <- exp$Tissue.Specificity
x4 <- exp$Duplication.Status
y <- dat.test$eigenvector[mt]

rm <- which(is.finite(log(y)) == FALSE | is.finite(log(x1)) == FALSE | is.finite(log(x2)) == FALSE)
mod_full <- lm(log(y[-rm]) ~ log(x1[-rm])*log(x2[-rm])*x3[-rm]*x4[-rm])
summary(mod_full)

mod_exp <- lm(log(y[-rm]) ~  x3[-rm] + log(x2[-rm]):x4[-rm])
summary(mod_exp)

plot(mod_exp$fitted.values, mod_exp$fitted.values + mod_exp$residuals)

quartz()
plot(log(x3[-rm]), log(y[-rm]))

#glmnet
dup <- model.matrix(~as.factor(x4[-rm]))[,-1] #dummy variables for WGD
reg.dat <- cbind(log(x1[-rm]), log(x2[-rm]), x3[-rm], dup)
colnames(reg.dat) <- c("Ka.KS", "Exp", "Tissue", "non_WGD", "Recent_WGD", "Singleton")

#standarizing (we don't really need to do this because glmnet standarizes, but it reports the un-standardized coefficients, so this makes it easier to compare)
reg.dat[,1:3] <- apply(reg.dat[,1:3], 2, function(x) return(x/sd(x)))

y.reg <- log(y[-rm])

alpha <- 1
# A tuning parameter for the lasso/elastic-net model-selection algorithm
# alpha = 1 corresponds to l1 penatly/lasso
# alpha = 0.5 corresponds to elastic net
# alpha = 0 gives ridge regression (l2 penalty), which will not do any model selection, merely shrinkage

#picking lambda (strenght of penalty by cross validation)
cv_linear_model <- cv.glmnet(x           = reg.dat,
                             y           =  y.reg,
                             family      = 'gaussian',
                             alpha       = alpha,
                             standardize = TRUE,
                             nfolds = 500)

best_lambda <- cv_linear_model$lambda.min

glmnet_mod <- glmnet(x = reg.dat, y = y.reg, family='gaussian', alpha=alpha, standardize=TRUE, lambda = best_lambda)

glmnet_mod$beta/sum(abs(glmnet_mod$beta))

#binning by centrality for a figure
IND <- rep("Peripheral", length(y.reg))
IND[which(y.reg > -30)] <- "Central"

boxplot(reg.dat[,"Tissue"] ~ IND, range = 0, col = "gray", ylab = "Degree of tissue specificity")
