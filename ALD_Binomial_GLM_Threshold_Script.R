#########################################################################
###
###
### MCMC Algorithm to fit binomial linear mixed model using quantile
### regression
###
###
#########################################################################

ALD_Binomial_GLM_MCMC=function(data,
                                tau=0.5,
                                inits,
                                parameters=c("alpha","beta"),
                                n.iter=5000,
                                checkpoint=1000,
                                out.loc="~/MCMC_Output.RData"
                                ){
    ##
    ##  Subroutines and Packages
    ##

    library(ald)

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

    p.ALD=function(q,mu=0,sigma=1,p=0.5){
        acum=ifelse(q<mu,
                    p*exp((1-p)*(q-mu)/sigma),
                    1-(1-p)*exp(-(p)*(q-mu)/sigma))
        return(acum)
    }

    q.ALD=function(prob,mu=0,sigma=1,p=0.5){
        q=ifelse(prob<p,
                 mu+((sigma*log(prob/p))/(1-p)),
                 mu-((sigma*log((1-prob)/(1-p)))/p))
        return(q)
    }

    r.truncALD=function(z,mu,tau,sigma=1){
        trunc=-mu/sigma
        g=p.ALD(q=trunc,mu=0,sigma=1,p=tau)
        min=ifelse(z==1,g,0)
        max=ifelse(z==1,1,g)
        u=min+(max-min)*runif(length(z))
        fn_val=ifelse(z==0&trunc<=0,
                      -(abs(trunc)-log(runif(length(z))))/(1-tau),
               ifelse(z==1&trunc>0,
                      trunc-log(runif(length(z)))/tau,
                      q.ALD(prob=u,mu=0,sigma=1,p=tau)
                      )
               )
        fn_val=mu+fn_val*sigma
        return(fn_val)
    }

    ##
    ## Priors
    ##

    beta.mean=priors$beta.prior[1]
    beta.var=priors$beta.prior[2]

    ##
    ## Dimensions
    ##

    n=length(data$y)
    p=length(inits$beta)

    ##
    ## Data
    ##

    y=data$y
    X=data$X

    ##
    ## Tuning parameters
    ##

    beta.tune=c(0.5,
                0.285,
                .2
                )
    beta.tune=rep(0.5,p)
    accept.beta=rep(0,p)
    t.tune=0.1
    accept.t=0

    ##
    ## Pre calculations
    ##

    beta=inits$beta
    mu=X%*%beta
    v=inits$v
    t=0

    ##
    ## Containers
    ##

    MCMC.Chains=vector('list',length(parameters))
    names(MCMC.Chains)=parameters
    if('beta'%in%parameters){
        MCMC.Chains$beta=matrix(,n.iter,p)
    }
    if('v'%in%parameters){
        MCMC.Chains$v=matrix(,n.iter,n)
    }
    if('mu'%in%parameters){
        MCMC.Chains$mu=matrix(,n.iter,n)
    }
    if('tuners'%in%parameters){
        MCMC.Chains$tuners=matrix(,n.iter,p)
    }
    if('t'%in%parameters){
        MCMC.Chains$t=matrix(,n.iter,1)
    }

    ##
    ## Gibbs loop
    ##

    for(k in 1:n.iter){

        ##
        ## Sample v
        ##

        v=r.truncALD(z=y,mu=mu,tau=tau,sigma=1)

        ##
        ## Sample beta
        ##

        for(i in 1:p){
            beta.tmp=rnorm(1,beta[i],beta.tune[i])
            beta.star=beta
            beta.star[i]=beta.tmp
            mu.star=X%*%beta.star
            mh1=sum(log(dALD.v(y=v,mu=mu.star,sigma=1,tau=tau))) +
                dnorm(beta.star[i],beta.mean,sqrt(beta.var),log=TRUE)
            mh2=sum(log(dALD.v(y=v,mu=mu,sigma=1,tau=tau))) +
                dnorm(beta[i],beta.mean,sqrt(beta.var),log=TRUE)
            mh=exp(mh1-mh2)
            if(min(mh,1)>runif(1)){
                beta=beta.star
                mu=mu.star
                accept.beta[i]=accept.beta[i]+1
            }
        }

        ##
        ## Sample t
        ##

        t.star=rnorm(1,t,t.tune)
        X.star=X
        X.star[,3]=X.star[,3]-t.star
        X.star[,3]=ifelse(X.star[,3]>t.star,X.star[,3]-t.star,0)
        mu.star=X.star%*%beta
        mh1=sum(log(dALD.v(y=v,mu=mu.star,sigma=1,tau=tau)))
        mh2=sum(log(dALD.v(y=v,mu=mu,sigma=1,tau=tau)))
        mh=exp(mh1-mh2)
        if(min(mh,1)>runif(1)){
            t=t.star
            X=X.star
            mu=mu.star
            accept.t=accept.t+1
        }

        ##
        ## Save samples
        ##

        if('beta'%in%parameters){
            MCMC.Chains$beta[k,]=beta
        }
        if('v'%in%parameters){
            MCMC.Chains$v[k,]=v
        }
        if('mu'%in%parameters){
            MCMC.Chains$mu[k,]=mu
        }
        if('t'%in%parameters){
            MCMC.Chains$t[k]=t
        }
        if('tuners'%in%parameters){
            MCMC.Chains$tuners[k,]=beta.tune
        }

        ##
        ## Checkpoint
        ##

        if(k%%checkpoint==0){

            ##
            ## Update tuning parameters
            ##

            beta.tune=ifelse(accept.beta/k>0.5,
                             beta.tune*1.1,
                      ifelse(accept.beta/k<0.3,
                             beta.tune*0.9,
                             beta.tune)
                      )
            if(accept.t/k>0.5){
                t.tune=t.tune*1.1
            }
            if(accept.t/k<0.3){
                t.tune=t.tune*0.9
            }

            ##
            ## Output results
            ##


        }
    }
    return(MCMC.Chains)
}
