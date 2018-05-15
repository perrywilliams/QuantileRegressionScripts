ALDOccupancyMCMC=function(data,
                          tau=0.5,
                          inits,
                          parameters=c("alpha","beta"),
                          n.iter=5000,
                          checkpoint=checkpoint,
                          out.loc){
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

    alpha.mean=priors$alpha.prior[1]
    alpha.var=priors$alpha.prior[2]
    beta.mean=priors$beta.prior[1]
    beta.var=priors$beta.prior[2]

    ##
    ## Tuning parameters
    ##

    alpha.tune= rep(0.01,length(alpha))
    beta.tune=rep(0.01,length(beta))
    accept.alpha=rep(0,length(alpha))
    accept.beta=rep(0,length(beta))

    ##
    ## Data and Dimensions
    ##

    y=as.matrix(data$y)
    n=dim(data$y)[1]
    m=length(inits$alpha)
    q=length(inits$beta)

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
    psi=1-p.ALD(-X%*%beta,mu=0,sigma=1,p=tau)
    p=1-p.ALD(-W%*%alpha,mu=0,sigma=1,p=tau)
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
    if('psi'%in%parameters){
        MCMC.Chains$psi=matrix(,n.iter,n)
    }
    if('p'%in%parameters){
        MCMC.Chains$p=matrix(,n.iter,n)
    }
    if('tuners'%in%parameters){
        MCMC.Chains$tuners=matrix(,n.iter,m+q)
    }

    ##
    ## Gibbs loop
    ##

    for(k in 1:n.iter){

        ##
        ## Sample z
        ##

        psi.temp=(pALD(psi,0,1,tau,TRUE)*
                  (1-pALD(p[,1],0,1,tau,TRUE))^J)/
            ((1-pALD(psi,0,1,tau,TRUE))+
             (pALD(psi,0,1,tau,TRUE)*
              (1-pALD(p[,1],0,1,tau,TRUE))^J)
            )
        z[y0]=rbinom(sum(y0),1,round(psi.temp[y0],3))
        z1=(z==1)
        z0=(z==0)

        ##
        ##  Sample alpha
        ##

        for(i in 1:length(alpha)){
            alpha.tmp=rnorm(1,alpha[i],alpha.tune[i])
            alpha.star=alpha
            alpha.star[i]=alpha.tmp
            p.star=matrix(1-p.ALD(-W%*%alpha.star,mu=0,sigma=1,p=tau),n,1)
            mh1=sum(dbinom(y[z1,],1,p.star[z1,],log=TRUE))+
                sum(dnorm(alpha.star[i],alpha.mean,sqrt(alpha.var),
                          log=TRUE))
            mh2=sum(dbinom(y[z1,],1,p[z1,],log=TRUE))+
                sum(dnorm(alpha[i],alpha.mean,sqrt(alpha.var),
                          log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>min(1,runif(1))){
                alpha=alpha.star
                p=p.star
                accept.alpha[i]=accept.alpha[i]+1
            }
        }

        ##
        ## Sample beta
        ##

        for(i in 1:length(beta)){
            beta.tmp=rnorm(1,beta[i],beta.tune[i])
            beta.star=beta
            beta.star[i]=beta.tmp
            psi.star=matrix(1-p.ALD(-X%*%beta.star,mu=0,sigma=1,p=tau),n,1)
            mh1=sum(dbinom(z,1,psi.star,log=TRUE))+
                sum(dnorm(beta.star[i],beta.mean,sqrt(beta.var),
                          log=TRUE))
            mh2=sum(dbinom(z,1,psi,log=TRUE))+
                sum(dnorm(beta[i],beta.mean,sqrt(beta.var),
                          log=TRUE))

            mh=exp(mh1-mh2)
            if(min(mh,1)>runif(1)){
                beta=beta.star
                psi=psi.star
                accept.beta[i]=accept.beta[i]+1
            }
        }


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
        if('psi'%in%parameters){
            MCMC.Chains$psi[k,]=psi
        }
        if('p'%in%parameters){
            MCMC.Chains$p[k,]=p[,1]
        }
        if('tuners'%in%parameters){
            MCMC.Chains$tuners[k,]=c(alpha.tune,beta.tune)
        }

        ##
        ## Checkpoint
        ##

        if(k%%checkpoint==0){

            ##
            ## Update tuning parameters
            ##

            alpha.tune=ifelse(accept.alpha/k>0.5,
                              alpha.tune*1.1,
                       ifelse(accept.alpha/k<0.3,
                              alpha.tune*0.9,
                              alpha.tune)
                       )

            beta.tune=ifelse(accept.beta/k>0.5,
                             beta.tune*1.1,
                      ifelse(accept.beta/k<0.3,
                             beta.tune*0.9,
                             beta.tune)
                      )

            ##
            ## Output results
            ##


        }
    }
    return(MCMC.Chains)
}
