## ---------------------------------------------------------------------------------------------------------
#===============================================================================
#===============================================================================
# Toy Example: Simulated data from an accelerated failure time model
#===============================================================================
#===============================================================================


################################################################################
# Required packages
################################################################################
rm(list=ls())
#library(devtools)
#install_github("davidrusi/mombf")
library(mombf)
#install_github("nyiuab/BhGLM")
library(BhGLM)
library(mvtnorm)
library(survival)
library(ggplot2)
library(glmnet)
library(BVSNLP)

################################################################################
# Data simulation
################################################################################

# Number of simulations (check n = 250 vs. n = 1000)
ns <- 250

# Simulation seed
set.seed(123)

# Number of active variables
pa <- 3

# number of spurious variables
ps <- 7

# Design matrix 
mu0 <- rep(0,pa+ps) # mean
Sigma0 <- diag(pa+ps) # Correlation matrix
Sigma0[lower.tri(Sigma0)] <- 0.5
Sigma0[upper.tri(Sigma0)] <- 0.5
X <- rmvnorm(n = ns, mean = mu0, sigma = Sigma0)
colnames(X) <- paste("Var",1:(pa+ps), sep = "")

# Scale all columns of X
X <- as.matrix(apply(X, 2, scale))

# Residuals
sd0 <- 0.5 # standard deviation
eps <- rnorm(n = ns, mean = 0, sd = sd0)

# True values of the regression coefficients
beta0 <- -0.25
betas <- c(0.25,0.5,1,rep(0,ps)) 

# Response variables
y0 <- as.vector(beta0 + X%*%betas + eps)

# Times to event
os <- exp(y0)

# Censoring times
cens <- rep(1, ns)

# Observed times
times <- ifelse(os  <= cens, os , cens)

# log(observed times)
ys <- log(times)

# Vital status indicators
status <-  ifelse(os  <= cens, 1 , 0)

# Data frame with all the data
df <- data.frame(times = times, ys = ys, status = status, X )

# Visualise the Kaplan-Meier estimator (of survival)
km <- survfit(Surv(times, status) ~ 1, data = df)
plot(km, xlab = "Years", ylab = "Survival", main = "Kaplan-Meier estimator",
     cex.axis = 1.5, cex.lab = 1.5)


## ---------------------------------------------------------------------------------------------------------
################################################################################
# Bayesian Variable Selection (BVS)
################################################################################

# Regression formula with linear effects
paste(colnames(X), collapse='+')
f1y <- as.formula(ys ~Var1+Var2+Var3+Var4+Var5+Var6+Var7+Var8+Var9+Var10)

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
priorCoefm= momprior(tau=0.192) # p-mom
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



## ---------------------------------------------------------------------------------------------------------
#####################################################################
# Cox PiMOM priors (sensitivity analysis)
#####################################################################

# Survival object with times
y2 <- Surv(time = df$times, event=df$status)

X2 <- X

bvsfit <- bvs(X = data.frame(X2), resp = as.matrix(cbind(times,status)), 
              family  = "survival", nlptype = "piMOM", mod_prior = "unif",
              niter = 100)

# Models with highest probability
head(bvsfit$max_models)

# Highest probabilities
head(exp(bvsfit$max_prob_vec)/sum(exp(bvsfit$max_prob_vec)))

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

