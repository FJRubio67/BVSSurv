## ---------------------------------------------------------------------------------------------------------
################################################################################
# Required packages
################################################################################
rm(list=ls())
#library(devtools)
#install_github("davidrusi/mombf")
library(mombf)
library(mvtnorm)
library(survival)
library(ggplot2)
library(glmnet)
library(penalized)
library(BhGLM)
library(BVSNLP)
library(eha)

################################################################################
# Data preparation
################################################################################

# nki70 data set
data(nki70)
head(nki70)
dim(nki70)
mean(nki70$event)

# Survival object
ys= Surv(log(nki70$time), event=nki70$event)
# Design matrix
X <- nki70[,-c(1:2)];
# Diam as binary
X$Diam <- as.numeric(X$Diam) - 1
# N as binary
X$N <- as.numeric(X$N)
# ER as binary
X$ER <- as.numeric(X$ER) - 1
# Grade as numeric
# Ideally, we would do Grade as factor
X$Grade <- as.numeric(X$Grade)

# Scale the remaining entries
X <- apply(X,2,scale)
#X[,-c(1:4)] <- apply(X[,-c(1:4)],2,scale)

df <- data.frame(ys = ys, as.data.frame(X), time = nki70$time, status = nki70$event)


# Visualise the Kaplan-Meier estimator (of survival)
km <- survfit(Surv(time, status) ~ 1, data = df)
plot(km, xlab = "Years", ylab = "Survival", main = "Kaplan-Meier estimator",
     cex.axis = 1.5, cex.lab = 1.5)



## ---------------------------------------------------------------------------------------------------------
################################################################################
# Bayesian Variable Selection (BVS)
################################################################################

# Regression formula with linear effects
f1y <- formula(paste('ys ~ ',paste(colnames(X), collapse='+'),sep=''))

#-------------------------------------------------------------------------------
# BVS: Prior 1 - Zellner's prior
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoef= zellnerprior(tau=nrow(X))
# Prior on the model space
priorDelta= modelbbprior(1,1)
# Selection step
ms1 <- modelSelection(f1y,data = df,
                      priorCoef=priorCoef,priorDelta=priorDelta,
                      enumerate=FALSE, method='Laplace', niter=10000, 
                      center = FALSE, scale = FALSE)

# Calculating model posterior probabilities
pp1 <- postProb(ms1)
# Top models
head(pp1)
# Marginal inclusion probabilities
mp1 <- ms1$margpp
names(mp1) <- c("Intercept",colnames(X))
# Visualising the marginal inclusion probabilities
df1<-data.frame(names=names(mp1),mp = as.vector(round(mp1,digits = 3)))
df1 <- df1[order(mp1, decreasing = TRUE),][1:10,]

ggplot(df1,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


# Fitting a CoxPH model using the survival R package
fit0 <- coxph(Surv(time, event) ~ PRC1 + KNTC2, data=nki70)
summary(fit0)

# AFT model using a lognormal baseline distribution
fit.aft <- aftreg(formula = Surv(time, event) ~ PRC1 + KNTC2, data=nki70,
                  dist = "lognormal")



#-------------------------------------------------------------------------------
# BVS: Prior 2 - MOM prior (Laplace approximation)
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoefm= momprior(taustd = 1) # p-mom
# Prior on the model space
priorDelta= modelbbprior(1,1)
# Selection step
ms2 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoefm,priorDelta=priorDelta,
                      enumerate=FALSE, method='Laplace', niter=10000,
                      center = FALSE, scale = FALSE)



# Calculating model posterior probabilities 
pp2 <- postProb(ms2)
# Top models
head(pp2)
# Marginal inclusion probabilities
mp2 <- ms2$margpp
names(mp2) <- c("Intercept",colnames(X))


# Visualising the marginal inclusion probabilities
df2<-data.frame(names=names(mp2),mp = as.vector(round(mp2,digits = 3)))
df2 <- df2[order(mp2, decreasing = TRUE),][1:10,]

ggplot(df2,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")

#-------------------------------------------------------------------------------
# BVS: Prior 3 - MOM prior (Approximate Laplace approximation)
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoefm= momprior(taustd = 1) # p-mom
# Prior on the model space
priorDelta= modelbbprior(1,1)
# Selection step
ms3 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoefm,priorDelta=priorDelta,
                      enumerate=FALSE, method='ALA', niter=10000,
                      center = FALSE, scale = FALSE)



# Calculating model posterior probabilities 
pp3 <- postProb(ms3)
# Top models
head(pp3)
# Marginal inclusion probabilities
mp3 <- ms3$margpp
names(mp3) <- c("Intercept",colnames(X))


# Visualising the marginal inclusion probabilities
df3 <-data.frame(names=names(mp3),mp = as.vector(round(mp3,digits = 3)))
df3 <- df3[order(mp3, decreasing = TRUE),][1:10,]

ggplot(df3,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")



## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Cox PiMOM priors (sensitivity analysis)
#####################################################################

# Survival object with times
y2 <- Surv(time = df$time, df$status)

bvsfit <- bvs(X = data.frame(X), resp = as.matrix(cbind(df$time,df$status)), 
              family  = "survival", nlptype = "piMOM", mod_prior = "unif",
              niter = 100)

# Models with highest probability
head(bvsfit$max_models)

# Highest probabilities
head(exp(bvsfit$max_prob_vec - max(bvsfit$max_prob_vec))/sum(exp(bvsfit$max_prob_vec - max(bvsfit$max_prob_vec) )))

# Highest posterior model
bvsfit$HPM

# Estimated coefficients
bvsfit$beta_hat

# Number of Visited Models
bvsfit$num_vis_models

# Visualising the marginal inclusion probabilities
dfi<-data.frame(names=names(mp1)[-1],mp = as.vector(round(bvsfit$inc_probs,digits = 3)))
dfi <- dfi[order(dfi$mp, decreasing = TRUE),][1:10,]

ggplot(dfi,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


# Fitting a CoxPH model using the survival R package
fit1 <- coxph(Surv(time, event) ~ PRC1 + KNTC2, data=nki70)
summary(fit1)




## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Spike-and-Slab LASSO (sensitivity analysis)
#####################################################################

# Selection step
spsl <- bmlasso(x = X, y =  y2, family = "cox")

# Posterior modes (only option in this package, no MCMC)
spsl$coef

# Posterior modes
plot.bh(coefs = spsl$coef, threshold = spsl$df, gap = 1, 
        main = "Spike-and-Slab LASSO", lty = 2) 

# Posterior Odds-Ratios for the modes
plot.bh(coefs = spsl$coef, threshold = spsl$df, gap = 1, 
        main = "Spike-and-Slab LASSO", lty = 2, OR = TRUE) 



## ---------------------------------------------------------------------------------------------------------
#############################################################################################
# Cox-LASSO (sensitivity analysis). Results as in https://doi.org/10.1002/sta4.607
#############################################################################################

# Selection step
cv.fit = try(cv.glmnet(x =  X, y = y2, family="cox", 
                       maxit=1e4, nfolds=10, alpha=1), silent=TRUE)
fit = try(glmnet(x = X, y=y2, family = "cox", maxit=1e4, alpha=1), silent=TRUE)

# active variables (lambda.min)
b.coxlasso = as.double(coef(fit, s=cv.fit$lambda.min))
which(b.coxlasso!=0)
colnames(X)[which(b.coxlasso!=0)]

# active variables (lambda.1se)
b2.coxlasso = as.double(coef(fit, s=cv.fit$lambda.1se))
which(b2.coxlasso!=0)
colnames(X)[which(b2.coxlasso!=0)]

# Fitting a CoxPH model using the survival R package
fit2 <- coxph(Surv(time, event) ~ N + QSCN6L1 + ZNF533 + IGFBP5.1 + PRC1, data=nki70)
summary(fit2)


