ProbitOccupancyMCMC=function(data,
                             inits,
                             parameters=c("alpha","beta"),
                             n.iter=5000,
                             checkpoint=checkpoint,
                             output.location){

    ##
    ##  Subroutines and Packages
    ##

    library(MASS)

    truncnormsamp <- function(mu,sig2,low,high,nsamp){
        flow=pnorm(low,mu,sqrt(sig2))
        fhigh=pnorm(high,mu,sqrt(sig2))
        u=runif(nsamp)
        tmp=flow+u*(fhigh-flow)
        x=qnorm(tmp,mu,sqrt(sig2))
        x
    }

    ##
    ## Dimensions
    ##

    n=dim(data$y)[1]
    m=length(inits$beta)
    J=data$J

    ##
    ## Pre calculations
    ##

    alpha=inits$alpha
    beta=inits$beta
    z=inits$z
    z1=(z==1)
    z0=(z==0)
    X=data$X
    W=data$W
    phi=X%*%beta
    p=W%*%alpha
    u=matrix(rnorm(n*max(J),p,1),n,max(J))
    v=matrix(rnorm(n,phi,1),n,1)
    y0=ifelse(apply(y,1,sum)==0,TRUE,FALSE)

    ##
    ## Containers
    ##

    MCMC.Chains=vector('list',length(parameters))
    names(MCMC.Chains)=parameters
    if('alpha'%in%parameters){
        MCMC.Chains$alpha=matrix(,n.iter,length(alpha))
    }
    if('beta'%in%parameters){
        MCMC.Chains$beta=matrix(,n.iter,length(beta))
    }
    if('z'%in%parameters){
        MCMC.Chains$z=matrix(,n.iter,n)
    }
    if('v'%in%parameters){
        MCMC.Chains$v=matrix(,n.iter,n)
    }
    if('phi'%in%parameters){
        MCMC.Chains$phi=matrix(,n.iter,n)
    }
    if('p'%in%parameters){
        MCMC.Chains$p=matrix(,n.iter,n)
    }

    ##
    ## Gibbs Loop
    ##

    for(k in 1:n.iter){

        ##
        ## Sample v
        ##

        v[z0]=truncnormsamp(phi[z0],1,-Inf,0,sum(z0))
        v[z1]=truncnormsamp(phi[z1],1,0,Inf,sum(z1))

        ##
        ## Sample beta
        ##

        var.tmp=solve(t(X)%*%X)
        mean.tmp=solve(t(X)%*%X)%*%t(X)%*%v
        beta=mvrnorm(1,mean.tmp,var.tmp)
        phi=X%*%beta

        ##
        ## Sample u
        ##

        y.mat=y[z1,]
        y.vec=c(y.mat)
        y.vec.0=(y.vec==0)
        p.vec=rep(p[z1],max(J))
        u.vec=rep(NA,length(y.vec))
        u.vec[y.vec.0]=truncnormsamp(p.vec[y.vec.0],1,-Inf,0,sum(y.vec.0))
        u.vec[!y.vec.0]=truncnormsamp(p.vec[!y.vec.0],1,0,Inf,sum(!y.vec.0))
        u[z1,]=u.vec

        ##
        ## Sample alpha
        ##

        W.tmp=W[z1,]
        for(i in 1:(max(J)-1)){
            W.tmp=rbind(W.tmp,W[z1,])
        }
        u.tmp=u[z1,]
        u.vec=c(u.tmp)
        var.tmp=solve(t(W.tmp)%*%W.tmp)
        mean.tmp=solve(t(W.tmp)%*%W.tmp)%*%t(W.tmp)%*%u.vec
        alpha=mvrnorm(1,mean.tmp,var.tmp)
        p=W%*%alpha

        ##
        ## Sample z
        ##

        phi.temp=(pnorm(phi,0,1)*(1-pnorm(p,0,1))^J)/
            (pnorm(phi,0,1)*(1-pnorm(p,0,1))^J+(1-pnorm(phi,0,1)))
        z[y0]=rbinom(sum(y0),1,round(phi.temp[y0],3))
        z1=(z==1)
        z0=(z==0)


        ##
        ## Save samples
        ##

        if('alpha'%in%parameters){
            MCMC.Chains$alpha[k,]=alpha
        }
        if('beta'%in%parameters){
            MCMC.Chains$beta[k,]=beta
        }
        if('z'%in%parameters){
            MCMC.Chains$z[k,]=z
        }
        if('v'%in%parameters){
            MCMC.Chains$v[k,]=v
        }
        if('phi'%in%parameters){
            MCMC.Chains$phi[k,]=phi
        }
        if('p'%in%parameters){
            MCMC.Chains$p[k,]=p
        }

        if(k%%checkpoint==0){
            save(MCMC.Chains,file=output.location)
        }
    }
}
