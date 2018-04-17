####################################################
###
### Probit Occupancy Model
###
### Created 01Mar2018
###
####################################################

rm(list=ls())

###
### Simulate occupancy data
###

set.seed(2017)
n=1000
m=4
alpha.truth=rnorm(m,0,3)
beta.truth=rnorm(m,0,3)
J=rep(3,n)
x1=rnorm(n)
x2=rnorm(n)
x3=x2^2
X=cbind(1,x1,x2,x3)
w1=rnorm(n)
w2=rnorm(n)
w3=x2^2
W=cbind(1,w1,w2,w3)
phi=X%*%beta.truth
v=rnorm(n,phi,1)
z=ifelse(v>0,1,0)
p=W%*%alpha.truth
u=matrix(rnorm(n*max(J),p,1),n,max(J))
y=matrix(0,n,max(J))
y.vec=c(y[z==1,])
u.vec=c(u[z==1,])
y.vec=ifelse(u.vec>0,1,0)
y[z==1,]=y.vec

data=list(y=y,
          X=X,
          W=W,
          J=J
          )

###
### MCMC Settings
###

n.iter=10000
checkpoint=1000

###
### Prior hyperparameters
###

## alpha
alpha.mean=0
alpha.var=10

## beta
beta.mean=0
beta.var=10

priors=list(alpha.prior=c(alpha.mean,alpha.var),
            beta.prior=c(beta.mean,beta.var)
            )


###
### Starting Values
###

alpha=alpha.truth
beta=beta.truth
z=matrix(ifelse(apply(y,1,sum)>0,1,0),n,1)

inits=list(alpha=alpha,
           beta=beta,
           z=z
           )

###
### Parameters to monitor
###

parameters=c("alpha",
             "beta",
             "p",
             "phi")

###
### Output location
###

output.location="~/Dropbox/QuantileRegression/Probit.RData"

###
### Fit model with 'ProbitOccupancyMCMC' function
###

source("~/Dropbox/GitHub/QuantileOccupancyModeling/ProbitOccupancyMCMC.R")

ProbitOccupancyMCMC(data=data,
                    inits=inits,
                    parameters=parameters,
                    n.iter=n.iter,
                    checkpoint=checkpoint,
                    output.location=output.location)

###
### Load and plot MCMC chains
###

load(output.location)
par(mfrow=c(4,2))
for(i in 1:4){
    plot(MCMC.Chains$alpha[,i],type='l')
    abline(h=alpha.truth[i],col=2)
}
for(i in 1:4){
    plot(MCMC.Chains$beta[,i],type='l')
        abline(h=beta.truth[i],col=2)
}

