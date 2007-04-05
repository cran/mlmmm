########################################################################
# PAN-EM example command file - works for Windows version of Splus4.5.
# Before running this example, make sure that the
# em.pan() and embd.pan() functions and object code are loaded into
# Splus; see the README file for details.
########################################################################
# The following dataset is Data Set 5.3, "Average Daily Gain", from the
# "SAS System for Mixed Models" by Little et. al.
# This dataset was collected to investigate average daily gains (AGD)
# of steers fed for 160 days. The treatments are four diets consisting
# of a base ration and three levels of a medicated feed additive added
# to the base ration. The objective of the experiment was to determine
# the optimal level of feed additive to maximize the average daily gain.
# The steers were housed in barns, the blocking factor, where each barn
# held four steers and the steers were individually fed.
########################################################################
# First we need to enter the response data into a matrix (or, as in
# this case, a vector) with one column for each variable.
y <- matrix(c(1.03,NA,
           1.54,477,
           1.82,444,
           NA,370,
           1.31,403,
           2.16,NA,
           2.13,450,
           NA,393,
           1.59,394,
           2.53,499,
           2.33,482,
           1.80,317,
           2.09,499,
           NA,411,
           2.21,391,
           2.82,396,
           1.66,371,
           2.30,418,
           2.65,486,
           2.18,NA,
           1.42,395,
           NA,325,
           1.58,316,
           1.49,311,
           1.41,414,
           1.65,313,
           NA,309,
           1.34,323,
           0.18,315,
           0.64,376,
           0.76,308,
           0.70,439),ncol=2,nrow=32,byrow=T)
########################################################################
# Next we create the vector subj(in this example subj is barn)
# whose length is equal to the number
# of rows in y. This vector indicates which rows of y belong to the
# first subject, which rows belong to the second subject, etc.
subj <- c(1,1,1,1,
       2,2,2,2,
       3,3,3,3,
       4,4,4,4,
       5,5,5,5,
       6,6,6,6,
       7,7,7,7,
       8,8,8,8)
#######################################################################
# For now, em.pan requires to specify a vector occasions, this
# requirement will be waived soon.
#
occ <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
#
########################################################################
# Now we must specify the model to be used for EM and/or imputation.
# Because the four treatments(????) are not clearly ordered,
# let's use a model that has
# an intercept and three dummy codes to allow the population means for
# the four occasions to be estimated freely. We will also allow the
# intercept to randomly vary by subject. For each subject i, the
# covariate matrices are then:
#
#                   1 1 0 0                  1
#                   1 0 1 0                  1
#           Xi =    1 0 0 1           Zi =   1
#                   1 0 0 0                  1
#
#
#
# When using pan() or pan.em(), these are combined into a single matrix called
# pred. The pred matrix has the same number of rows as y, the number of
# subject-occasions. Each column of Xi and Zi must be represented in
# pred. Because Zi is merely the first column of Xi, we do not need to
# enter that column twice. So pred is simply the matrix Xi, stacked
# upon itself nine times.
#
pred <- cbind(int=rep(1.,32),
   dummy1=rep(c(1.,0.,0.,0.),8),
   dummy2=rep(c(0.,1.,0.,0.),8),
   dummy3=rep(c(0.,0.,1.,0.),8))
#
# Now we must tell pan that all six columns of pred are to be used in
# Xi, but only the first column of pred appears in Zi.
#
xcol <- 1:4
zcol <- 1
########################################################################
# The model specification is now complete.
# Now we are ready to run em.pan(). Let's assume that the em.pan function
# and the object code have already been loaded into Splus.
#
res <- em.pan(y,subj,pred,xcol,zcol,maxits=200,eps=0.0001)
#
# Now, we have ML estimates of the model parameters beta, Psi and Sigma.
# These estimates can be used as prior guesttimates for variance
# components when the missing values are imputed via pan() function.
# Recall that the dimension of Sigma is (r x r) where r
# is the number of response variables (in this case, r=1). The prior
# distribution for Sigma is inverted Wishart with hyperparameters a
# (scalar) and Binv (r x r), where a is the imaginary degrees of freedom
# and Binv/a is the prior guesstimate of Sigma. The value of a must be
# greater than or equal to r. The "least informative" prior possible
# would have a=r, so here we will take a=2. As a prior guesstimate of
# Sigma we will use the ML estimate of Sigma from "res", res$sigma, therefore
# Binv=2*res$sigma
#
# By similar reasoning we choose the prior distribution for Psi. The
# dimension of Psi is (r*q x r*q) where q is the number of random
# effects in the model (i.e. the length of zcol, which in this case is
# one). The hyperparameters for Psi are c and Dinv, where c is the
# imaginary degrees of freedom (which must be greater than or equal to
# r*q) and res$psi is the prior guesstimate of Psi. We will take d=2
# and Cinv=2*res$psi.
#
# The prior is specified as a list with four components named a, Binv,
# c, and Dinv, respectively.
#
prior <- list(a=2,Binv=2*res$sigma,c=2,Dinv=2*res$psi)
########################################################################
# Now we are ready to run pan(). Let's assume that the pan function
# and the object code have already been loaded into Splus. First we
# do a preliminary run of 500 iterations.
#
result <- pan(y,subj,pred,xcol,zcol,prior,seed=13579,iter=1000)
#
# Check the convergence behavior by making time-series plots and acfs
# for the model parameters. Variances will be plotted on a log
# scale. We'll assume that a graphics device has already been opened.
#
plot(1:1000,log(result$sigma[1,1,]),type="l")
acf(log(result$sigma[1,1,]))
plot(1:1000,log(result$psi[1,1,]),type="l")
acf(log(result$psi[1,1,]))
par(mfrow=c(4,2))
for(i in 1:4) plot(1:1000,result$beta[i,1,],type="l")
for(i in 1:4) acf(result$beta[i,1,])
#
# This example appears to converge very rapidly; the only appreciable
# autocorrelations are found in Beta, and even those die down by lag
# 20. With a sample this small we can afford to be cautious, so let's
# impute the missing data m=10 times taking 100 steps between
# imputations. We'll use the current simulated value of y as the first
# imputation, then restart the chain where we left off to produce
# the second thru the tenth.
#
#
y1 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9565,iter=100,start=result$last)
y2 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=6047,iter=100,start=result$last)
y3 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=3955,iter=100,start=result$last)
y4 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=4761,iter=100,start=result$last)
y5 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9188,iter=100,start=result$last)
y6 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9029,iter=100,start=result$last)
y7 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=4343,iter=100,start=result$last)
y8 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=2372,iter=100,start=result$last)
y9 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=7081,iter=100,start=result$last)
y10 <- result$y
########################################################################
