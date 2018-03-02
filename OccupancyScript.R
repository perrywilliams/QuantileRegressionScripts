## ##########################################################################
## ###
## ### Quantile regression for occupancy surveys
## ###
## ### Contents:
## ### 1) Simulated data analysis with probit occupancy
## ### 2) Simulated data analysis with ALD occupancy
## ### 3) Willow tit analysis code
## ###
## ##########################################################################


## ####################################################
## ### 1) Simulated data analysis with probit occupancy
## ####################################################

## rm(list=ls())

## ###
## ### Simulate occupancy data
## ###

## n=100
## m=4
## p.truth=rbeta(1,1,1)
## beta.truth=rnorm(m,0,3)
## J=sample(c(1,2,3),n,replace=TRUE,prob=c(0.1,0.25,0.65))
## x1=rnorm(n)
## x2=rnorm(n)
## x3=x2^2
## X=cbind(1,x1,x2,x3)
## m=dim(X)[2]
## phi=rnorm(n,X%*%beta.truth,1)
## z=ifelse(phi<0,0,1)
## y=ifelse(z==1,rbinom(n,J,p.truth),0)

## ###
## ### MCMC Algorithm
## ###

## n.iter=5000

## ###
## ###  Subroutines
## ###

## truncnormsamp <- function(mu,sig2,low,high,nsamp){
##     flow=pnorm(low,mu,sqrt(sig2))
##     fhigh=pnorm(high,mu,sqrt(sig2))
##     u=runif(nsamp)
##     tmp=flow+u*(fhigh-flow)
##     x=qnorm(tmp,mu,sqrt(sig2))
##     x
## }

## ###
## ### Prior hyperparameters
## ###

## ## p
## p.alpha=1
## p.beta=1

## ## beta
## beta.var=10^2
## Sig0=beta.var*diag(m)
## beta.mu=matrix(0,m,1)


## ###
## ### Starting Values
## ###

## beta=beta.truth
## z=ifelse(y>0,1,0)
## z1=(z==1)
## z0=(z==0)

## ###
## ### Containers
## ###

## p.save=numeric(n.iter)
## beta.save=matrix(,n.iter,m)

## ###
## ### Fit model with MCMC
## ###

## for(k in 1:n.iter){

##     ##
##     ##  Sample p
##     ##

##     alpha.tmp=p.alpha+sum(y[z==1])
##     beta.tmp=p.beta+sum((J-y)[z==1])
##     p=rbeta(1,alpha.tmp,beta.tmp)

##     ##
##     ## Sample z
##     ##

##     phi.temp=(pnorm(phi)*(1-p)^J)/(pnorm(phi)*(1-p)^J+(1-pnorm(phi)))
##     z[y==0]=rbinom(sum(y==0),1,round(phi.temp[y==0],3))
##     z1=(z==1)
##     z0=(z==0)


##     ##
##     ## Sample phi
##     ##

##     phi=matrix(0,n,1)
##     phi1=truncnormsamp((X%*%beta)[z1],1,0,Inf,sum(z))
##     phi0=truncnormsamp((X%*%beta)[z0],1,-Inf,0,sum(1-z))
##     phi[z1,]=phi1
##     phi[z0,]=phi0
##     ## z=ifelse(phi<0,0,1)

##     ##
##     ## Sample beta
##     ##

##     var.tmp=solve(t(X)%*%X + solve(Sig0))
##     mean.tmp=var.tmp%*%(t(X)%*%phi + solve(Sig0)%*%beta.mu)
##     beta=mean.tmp+t(chol(var.tmp))%*%matrix(rnorm(m),m,1)

##     ##
##     ## Save samples
##     ##

##     p.save[k]=p
##     beta.save[k,]=beta

##     out=list(p.save,
##              beta.save)

## }

## par(mfrow=c(3,2))
## plot(out[[1]],type='l')
## abline(h=p.truth,col=2)

## for(i in 1:m){
##     plot(out[[2]][,i],type='l')
##     abline(h=beta.truth[i],col=2)
## }

####################################################
####################################################
###
### 2) Simulated data analysis with ALD occupancy
###
####################################################
####################################################

rm(list=ls())
## install.packages("ald")
library(ald)

###
###  Subroutines
###

dALD.v=function(y,mu,sigma,tau){
    densi=ifelse(y<mu,
    (tau*(1-tau)/sigma)*exp((1-tau)*(y-mu)/sigma),
    (tau*(1-tau)/sigma)*exp(-(tau)*(y-mu)/sigma)
    )
    return(densi)
}

dALD.m=function(y,mu,sigma,tau){
    densi=matrix(ifelse(test=c(y)<c(mu),
                        (tau*(1-tau)/sigma)*exp((1-tau)*(c(y)-c(mu))/sigma),
                        (tau*(1-tau)/sigma)*exp(-(tau)*(c(y)-c(mu))/sigma)),
                 dim(y)[1],dim(y)[2])
    return(densi)
}


rtruncALD=function(z,mu,tau,sigma=1){
    trunc=-mu/sigma
    if((z==0)&(trunc<=0)){
        u=runif(1)
        fn_val=-(abs(trunc)-log(u))/(1-tau)
    }else if((z==1)&(trunc>0)){
        u=runif(1)
        fn_val=trunc-log(u)/tau
    }else{
        g=pALD(q=trunc,mu=0,sigma=1,p=tau)
        if(z==1){
            min=g;max=1
        }else{
            min=0;max=g
        }
        u=min+(max-min)*runif(1)
        fn_val=qALD(prob=u,mu=0,sigma=1,p=tau)
    }
    fn_val=mu+fn_val*sigma
    return(fn_val)
}

###
### Simulate occupancy data
###

n=500 # no. sites
m=4   # no. beta parameters
alpha.truth=c(2,1)
beta.truth=c(2,1,3,4)
sigma=1
J=rep(3,n) ## no. visits per site
x1=rnorm(n)
x2=rnorm(n)
x3=rnorm(n)
X=cbind(1,x1,x2,x3)
w1=rnorm(n)
W=cbind(1,w1)
phi.truth=X%*%beta.truth
p.truth=W%*%alpha.truth
v.truth=sapply(1:n,function(i)rALD(1,mu=phi.truth[i],sigma=1,p=0.5))
z.truth=ifelse(v.truth<0,0,1)
z1=(z.truth==1)
z0=(z.truth==0)
u.truth=matrix(NA,n,max(J))
for(i in 1:3){
    u.truth[,i]=sapply(1:n,function(i)rALD(1,mu=p.truth[i],sigma=1,p=0.5))
}
y=matrix(NA,dim(u.truth)[1],dim(u.truth)[2])
for(i in 1:max(J)){
     y[z1,i]=ifelse(u.truth[z1,i]<0,0,1)
}
y[z0,]=0
zero.y=ifelse(apply(y,1,sum)>0,FALSE,TRUE)
truth=list()
truth[[1]]=alpha.truth
truth[[2]]=beta.truth
truth[[3]]=phi.truth
truth[[4]]=v.truth
truth[[5]]=u.truth
truth[[6]]=z.truth
truth[[7]]=y
names(truth)=c("alpha","beta","phi","v","u","z","y")

###
### MCMC Algorithm
###

n.iter=20000
checkpoint=100

###
### Prior Parameters
###

## alpha
alpha.var=10^2
Sig0.alpha=alpha.var*diag(1)
alpha.mu=matrix(0,1,1)

## beta
beta.var=10^2
Sig0.beta=beta.var*diag(m)
beta.mu=matrix(0,m,1)

###
### Starting Values
###

beta=beta.truth
phi=X%*%beta
alpha=alpha.truth
p=matrix(W%*%alpha,n,3)
u=u.truth
v=v.truth
z=z.truth
z1=(z==1)
z0=(z==0)
tau=0.5
accept.alpha=0
accept.beta=rep(0,m)

###
### Containers
###

alpha.save=matrix(,n.iter,2)
alpha.tune.save=numeric(n.iter)
beta.save=matrix(,n.iter,m)
beta.tune.save=matrix(NA,n.iter,m)

###
### Tuning parameters
###

alpha.tune= 0.4237573
beta.tune=c(0.1899601,0.2321735,0.4457457,0.1718543)

###
### Fit model with MCMC
###

for(k in 1:n.iter){

    ##
    ## Sample z
    ##

    phi.temp=(pALD(phi,0,sigma,tau,TRUE)*
              (1-pALD(p[,1],0,sigma,tau,TRUE))^J)/
        (
            (1-pALD(phi,0,sigma,tau,TRUE))+
            (pALD(phi,0,sigma,tau,TRUE)*(1-pALD(p[,1],0,sigma,tau,TRUE))^J)
        )
    z[zero.y]=rbinom(sum(zero.y),1,round(phi.temp[zero.y],3))
    z1=(z==1)
    z0=(z==0)

    ##
    ## Sample v
    ##

    v=sapply(1:n,function(i)rtruncALD(
                                z=z[i],
                                mu=phi[i],
                                tau=tau,
                                sigma=sigma))

    ##
    ## Sample u
    ##

    y.mat=y[z1,]
    y.vec=c(y.mat)
    p.mat=p[z1,]
    p.vec=c(p.mat)
    u.vec=rep(NA,length(y.vec))
    u.vec=sapply(1:length(u.vec),
                          function(i)rtruncALD(z=y.vec[i],
                                               mu=p.vec[i],
                                               tau=tau,
                                               sigma=sigma))
    u[z1,]=u.vec

    ##
    ##  Sample alpha
    ##

    alpha.star=rnorm(length(alpha),alpha,alpha.tune)
    p.star=matrix(W%*%alpha.star,n,3)
    mh1=sum(log(dALD.m(y=u[z1,],mu=p.star[z1,],sigma=sigma,tau=tau)))
    mh2=sum(log(dALD.m(y=u[z1,],mu=p[z1,],sigma=sigma,tau=tau)))
    mh=exp(mh1-mh2)
    if(mh>min(1,runif(1))){
        alpha=alpha.star
        p=p.star
        accept.alpha=accept.alpha+1
    }

    ##
    ## Sample beta0
    ##

    beta0.star=rnorm(1,beta[1],beta.tune[1])
    beta.star=c(beta0.star,beta[2:m])
    phi.star=X%*%beta.star
    mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
    mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
    mh=exp(mh1-mh2)
    if(min(mh,1)>runif(1)){
        beta=beta.star
        phi=phi.star
        accept.beta[1]=accept.beta[1]+1
    }

    ##
    ## Sample beta1
    ##

    beta1.star=rnorm(1,beta[2],beta.tune[2])
    beta.star=c(beta[1],beta1.star,beta[3:m])
    phi.star=X%*%beta.star
    mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
    mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
    mh=exp(mh1-mh2)
    if(min(mh,1)>runif(1)){
        beta=beta.star
        phi=phi.star
        accept.beta[2]=accept.beta[2]+1
    }

    ##
    ## Sample beta2
    ##

    beta2.star=rnorm(1,beta[3],beta.tune[3])
    beta.star=c(beta[1:2],beta2.star,beta[4:m])
    phi.star=X%*%beta.star
    mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
    mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
    mh=exp(mh1-mh2)
    if(min(mh,1)>runif(1)){
        beta=beta.star
        phi=phi.star
        accept.beta[3]=accept.beta[3]+1
    }

    ##
    ## Sample beta3
    ##

    beta3.star=rnorm(1,beta[4],beta.tune[4])
    beta.star=c(beta[1:3],beta3.star)
    phi.star=X%*%beta.star
    mh1=sum(log(dALD.v(y=v,mu=phi.star,sigma=sigma,tau=tau)))
    mh2=sum(log(dALD.v(y=v,mu=phi,sigma=sigma,tau=tau)))
    mh=exp(mh1-mh2)
    if(min(mh,1)>runif(1)){
        beta=beta.star
        phi=phi.star
        accept.beta[4]=accept.beta[4]+1
    }

    ##
    ## Save samples
    ##

    alpha.save[k,]=alpha
    alpha.tune.save[k]=alpha.tune
    beta.save[k,]=beta
    beta.tune.save[k,]=beta.tune

    ##
    ## Checkpoint
    ##

    if(k%%checkpoint==0){

        ##
        ## Update tuning parameters
        ##

        if(accept.alpha/k<0.3){
            alpha.tune=alpha.tune*0.9
        }
        if(accept.alpha/k>0.5){
            alpha.tune=alpha.tune*1.1
        }

        beta.tune=ifelse(accept.beta/k>0.5,
                         beta.tune*1.1,
                  ifelse(accept.beta/k<0.3,
                         beta.tune*0.9,
                         beta.tune)
                  )

        ##
        ## Export results
        ##

        out=list(alpha.save,
                 alpha.tune.save,
                 beta.save,
                 beta.tune.save,
                 truth)

        save(out,file=paste("~/Dropbox/QuantileRegression/",
                            "ALDData.ALDModel_6.RData",sep="")
             )
    }
}


###########################################################################
###
### Results Files:
###        ProbitData.ALDModel_1.RData (26 Feb 2018):
###               sample betas,alpha, and v only (didn't work).
###        ProbitData.ALDModel_2.RData (26 Feb 2018):
###               sample beta0 and alpha only (worked).
###        ProbitData.ALDModel_3.RData (26 Feb 2018):
###               sample beta0 alpha and v only (didn't work).
###        ProbitData.ALDModel_4.RData (26 Feb 2018):
###               sample beta0,beta1,alpha,v only (didn't work).
###               Note: may not have been working becuase different
###               model used to generate the data.
###               Now: switching to ALD generating model.
###        ALDData.ALDModel_1.RData (26 Feb 2018):
###               sample betas and v only. (worked)
###        ALDData.ALDModel_2.RData (26 Feb 2018):
###               sample betas and v only. (worked)
###        ALDData.ALDModel_3.RData (26 Feb 2018):
###               sample betas, v, z only. (worked)
###        ALDData.ALDModel_4.RData (27 Feb 2018):
###               sample alpha, v, z, and u (worked).
###        ALDData.ALDModel_5.RData (27 Feb 2018):
###               sample alpha, beta, v, z, and u (worked)
###        ALDData.ALDModel_6.RData (27 Feb 2018):
###               like 5, but two alphas
###
###########################################################################

rm(list=ls())
## Load and plot results


load(paste("~/Dropbox/QuantileRegression/",
           "ALDData.ALDModel_6.RData",sep="")
     )

(status=sum(!is.na(out[[1]][,1])))
burn=1

###
### Plot results
###

## alpha
par(mfrow=c(2,2))
plot(out[[1]][burn:status,1],type='l')
abline(h=out[[5]]$alpha[1],col=2)
plot(out[[2]][burn:status],type='l')
plot(out[[1]][burn:status,2],type='l')
abline(h=out[[5]]$alpha[2],col=2)
plot(out[[2]][burn:status],type='l')

## beta
par(mfrow=c(4,2))
beta.truth=c(2,1,3,4)
for(i in 1:4){
    plot(out[[3]][burn:status,i],type='l')
    abline(h=beta.truth[i],col=2)
    plot(out[[4]][burn:status,i],type='l')
}
