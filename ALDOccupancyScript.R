####################################################
###
### ALD Occupancy Model
###
### Created 01Mar2018
###
####################################################

rm(list=ls())
library(ald)
library(RCurl)

###
### Quantile
###

tau=0.5

###
### Simulate occupancy data
###

n=500
m=4
alpha.truth=rnorm(4)
beta.truth=rnorm(4)
sigma=1
J=rep(3,n)
x1=rnorm(n)
x2=rnorm(n)
x3=rnorm(n)
X=cbind(1,x1,x2,x3)
w1=rnorm(n)
w2=rnorm(n)
w3=x2^2
W=cbind(1,w1,w2,w3)
phi.truth=X%*%beta.truth
v.truth=sapply(1:n,function(i)rALD(1,mu=phi.truth[i],sigma=1,p=tau))
z.truth=ifelse(v.truth<0,0,1)
p.truth=W%*%alpha.truth
u.truth=matrix(NA,n,max(J))
for(i in 1:3){
    u.truth[,i]=sapply(1:n,function(i)rALD(1,mu=p.truth[i],sigma=1,p=tau))
}
y=matrix(0,dim(u.truth)[1],dim(u.truth)[2])
for(i in 1:max(J)){
     y[z.truth==1,i]=ifelse(u.truth[z.truth==1,i]<0,0,1)
}
y.vec=c(y[z.truth==1,])
u.vec=c(u.truth[z.truth==1,])
y.vec=ifelse(u.vec>0,1,0)
y[z.truth==1,]=y.vec

###
### Bundle data for MCMC sampler
###

data=list(y=y,
          X=X,
          W=W,
          J=J
          )

###
### MCMC Settings
###

n.iter=50000
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
             "beta")

###
### Output location
###

output.location="~/Output.RData"

###
### Fit model with 'ALDOccupancyMCMC' function
###

script=getURL(
    paste("https://raw.githubusercontent.com/",
          "perrywilliams/QuantileRegressionScripts/",
          "master/ALDOccupancyMCMC.R",sep=""),
    ssl.verifypeer=FALSE)
eval(parse(text = script))

sys.time=Sys.time()
ALDOccupancyMCMC(data=data,
                 tau=0.5,
                 inits=inits,
                 parameters=parameters,
                 n.iter=n.iter,
                 checkpoint=checkpoint,
                 output.location=output.location)
Sys.time()-sys.time

###
### Load and plot MCMC chains
###

load(output.location)
par(mfrow=c(4,2))
for(i in 1:4){
    plot(MCMC.Chains$alpha[,i],type='l',
         ylab="")
    abline(h=alpha.truth[i],col=2)
}
for(i in 1:4){
    plot(MCMC.Chains$beta[,i],type='l',
         ylab="")
        abline(h=beta.truth[i],col=2)
}
