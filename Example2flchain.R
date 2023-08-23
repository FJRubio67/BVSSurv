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
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biobase")
library(Biobase)
library(BhGLM)
library(BVSNLP)
library(eha)

################################################################################
# Data preparation
################################################################################

# flchain data set
data(flchain)
head(flchain)

# Variables of interest
colnames(flchain)[c(2,4,5,7,8,9,10)]

# Complete cases of the variables of interest
data.c <- flchain[,c(2,4,5,7,8,9,10)]
data.c <- data.c[complete.cases(data.c),]
dim(data.c)
head(data.c)

# Assigning half-day survival to the zero-survivors
data.c$futime <- ifelse(data.c$futime==0,0.5,data.c$futime)
ys <- Surv(log(data.c$futime/365.25), data.c$death)
data.c$kappa <- scale(data.c$kappa) # scaled kappa
data.c$creatinine <- scale(data.c$creatinine) # scaled creatinine 
data.c$lambda <- scale(data.c$lambda) # scaled lambda
data.c$mgus <- scale(data.c$mgus) # scaled mgus 
data.c$sex <- as.numeric(data.c$sex)-1 # sex 0-1
data.c$sex <- scale(data.c$sex) # scaled sex 
X <- data.c[,-c(6,7)];  # Design matrix
dim(X)

# Data frame with all the data
df <- data.c

# Visualise the Kaplan-Meier estimator (of survival)
km <- survfit(Surv(futime/365.25, death) ~ 1, data = df)
plot(km, xlab = "Years", ylab = "Survival", main = "Kaplan-Meier estimator",
     cex.axis = 1.5, cex.lab = 1.5)


## ---------------------------------------------------------------------------------------------------------
################################################################################
# Bayesian Variable Selection (BVS)
################################################################################

# Regression formula with linear effects
paste(colnames(X), collapse='+')
f1y <- as.formula(ys ~sex+kappa+lambda+creatinine+mgus)

#-------------------------------------------------------------------------------
# BVS: Prior 1 - Zellner's prior
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoef= zellnerprior(tau=nrow(X))
# Prior on the model space
priorDelta= modelbbprior(1,1)
# Selection step
ms1 <- modelSelection(f1y,data=df,
                       priorCoef=priorCoef,priorDelta=priorDelta,
                       enumerate=TRUE, method='Laplace', 
                       center = FALSE, scale = FALSE)

# Calculating model posterior probabilities
pp1 <- postProb(ms1)
# Top models
head(pp1)
# Marginal inclusion probabilities
mp1 <- ms1$margpp

# Visualising the marginal inclusion probabilities
df1<-data.frame(names=names(mp1),mp = as.vector(round(mp1,digits = 3)))
df1

ggplot(df1,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


#-------------------------------------------------------------------------------
# BVS: Prior 2 - MOM prior (Laplace)
#-------------------------------------------------------------------------------
# Priors for the coefficients
priorCoefm= momprior(tau = 0.192) # p-mom
# Prior on the model space
priorDelta= modelbbprior(1,1)
# Selection step
ms2 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoefm,priorDelta=priorDelta,
                      enumerate=TRUE, method='Laplace', 
                      center = FALSE, scale = FALSE)



# Calculating model posterior probabilities 
pp2 <- postProb(ms2)
# Top models
head(pp2)
# Marginal inclusion probabilities
mp2 <- ms2$margpp


# Visualising the marginal inclusion probabilities
df2<-data.frame(names=names(mp2),mp = as.vector(round(mp2,digits = 3)))
df2

ggplot(df2,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


#-------------------------------------------------------------------------------
# BVS: Prior 2 - MOM prior (ALA)
#-------------------------------------------------------------------------------

# Selection step
ms3 <- modelSelection(f1y,data=df,
                      priorCoef=priorCoefm,priorDelta=priorDelta,
                      enumerate=TRUE, method='ALA', 
                      center = FALSE, scale = FALSE)



# Calculating model posterior probabilities 
pp3 <- postProb(ms3)
# Top models
head(pp3)
# Marginal inclusion probabilities
mp3 <- ms3$margpp


# Visualising the marginal inclusion probabilities
df3<-data.frame(names=names(mp3),mp = as.vector(round(mp3,digits = 3)))
df3

ggplot(df3,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")

# AFT model using a lognormal baseline distribution
fit.aft <- aftreg(formula = Surv(futime, death) ~ kappa + lambda, data=df,
                  dist = "lognormal")

summary(fit.aft)

## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Cox PiMOM priors (sensitivity analysis)
#####################################################################

# Survival object with times
y2 <- Surv(time = data.c$futime/365.25, data.c$death)

X2 <- X

bvsfit <- bvs(X = data.frame(X2), resp = as.matrix(cbind(data.c$futime/365.2,data.c$death)), 
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
dfi

ggplot(dfi,aes(x =  reorder(names, -mp), y = mp,label=mp))+
  geom_bar(stat="identity")+geom_text(size=2.5,hjust=0)+coord_flip()+
  xlab("") + ylab("Marginal inclusion probability")


## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Spike-and-Slab LASSO (sensitivity analysis)
#####################################################################

# Selection step
spsl <- bmlasso(x = X2, y =  y2, family = "cox")

# Posterior modes (only option in this package, no MCMC)
spsl$coef

# Posterior modes
plot.bh(coefs = spsl$coef, threshold = spsl$df, gap = 1, 
        main = "Spike-and-Slab LASSO", lty = 2) 

# Posterior Odds-Ratios for the modes
plot.bh(coefs = spsl$coef, threshold = spsl$df, gap = 1, 
        main = "Spike-and-Slab LASSO", lty = 2, OR = TRUE) 


## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Cox-LASSO (sensitivity analysis)
#####################################################################

# Survival object with times
y2 <- Surv(time = df$futime/365.25, event=df$death)

X2 <- as.matrix(X)

# Selection step
cv.fit = try(cv.glmnet(x =  X2, y = y2, family="cox", 
                       maxit=10000, nfolds=10, alpha=1), silent=TRUE)
fit = try(glmnet(x = X2, y=y2, family = "cox", maxit=10000, alpha=1), silent=TRUE)

# active variables (lambda.min)
b.coxlasso = as.double(coef(fit, s=cv.fit$lambda.min))
which(b.coxlasso!=0)
colnames(X)[which(b.coxlasso!=0)]

# active variables (lambda.1se)
b2.coxlasso = as.double(coef(fit, s=cv.fit$lambda.1se))
which(b2.coxlasso!=0)
colnames(X)[which(b2.coxlasso!=0)]

# Fitting a CoxPH model using the survival R package
fit2 <- coxph(Surv(futime, death) ~ kappa + lambda, data=df)
summary(fit2)
