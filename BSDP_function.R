#### library import
library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma);library(bsplinePsd)


####
#function
####

### define the function of DPW: stick-breaking  
weights2=function(K,L,V,Z){
  p = rep(NA, L)
  p[1] = V[1]
  for (i in 2:L) {
    p[i] = prod(1 - V[1:(i - 1)]) * V[i]
  }
  p = c(1 - sum(p), p)  
  
  Zbound = matrix(0,K,2)
  for (i in 1:K) {
    Zbound[i,] = c((i-1)/K,i/K)
  }
  Zg = rep(0,L+1)
  # weights of each beta density
  G = rep(0, K) 
  for (j in 1:L+1){
    for (i in 1:nrow(Zbound)) {
      if ((Z[j]>Zbound[i,1]) & (Z[j]<=Zbound[i,2])){
        G[i]=G[i]+p[j]
        Zg[j]=i
        break
      }
    }
  }
  out=list(G=G,Zg=Zg,P=p)
  return(out)
}


### the function of posterior computation
#suppose both K,M in G and Bspline are the same
# lambda's prior is Gamma(a,b)
# sigma's prior is IGamma(nu0,s0)
# for Y's range use mu_f, baseb is the base
# with burnin 2000 iterations
# dk is the MH parameter of K

gibbs_bspline_est_linkfunction= function(
    X,Y,S=10000,L=max(20,n^(1/3)),Kmax=50,M=1,a=0.01,b=0.01,K=10,d=3,nu0=500,s0=0.05,lambda=10,dk=4,sigma=0.5,baseb=1.15,burnin=2000,logstate=TRUE
){
    n=nrow(X)
    p=ncol(X)

    ###
    # The initialization
    ###

    # V Z for the weights
    G=rep(NA,K)
    V = rep(1/2,L)
    Z = runif(L + 1)
    DPW=weights2(K,L,V,Z)
    G=DPW$G

    # U C for the knots
    U = rep(1/2,L)
    C = runif(L + 1)
    DPK=weights2(K-d,L,U,C)
    knot_diff=c(0,DPK$G)
    in_knots=rep(0,(K-d+1))
    in_knots=cumsum(knot_diff)

    # Initialise r
    r=rbinom(p,1,1/2)
    if(all(r==0)){r[sample(1:p,1)]=1}

    #Initialise beta
    mcmc_r_beta_ratio=rexp(1,3)
    mcmc_beta_r_rho=0.2
    beta=rep(0,p)
    beta[which(r==1)]=rnorm(length(which(r==1)),0,1) #unit norm
    beta=beta/c(sqrt(crossprod(beta)));
    if(beta[which(r==1)[1]]<0){beta=-beta} #beta1>0

    #Initialise kc
    kc=runif(n,0,3) 
    
    #if log
    if(logstate){
      mu_f=function(x,baseb){
          return(logb(x,baseb))
        } 
    }else{
      mu_f=function(x,baseb){
            return(x)
        }
    }
    
    ###
    # Parameter space
    ###
    KK=Lambda=bspline_dens=Raa=Sigma=NULL
    R=Beta=matrix(nr=S,nc=p)
    VV=matrix(nr=S,nc=L)
    ZZ=matrix(nr=S,nc=L+1)
    UU=matrix(nr=S,nc=L)
    CC=matrix(nr=S,nc=L+1)
    KC=matrix(nr=S,nc=n)
    GG=matrix(nr=S,nc=Kmax)
    GG_knots=matrix(nr=S,nc=Kmax)
    Rho=R_beta_ratio=NULL

    ###
    # target parameters updating process
    ###
    for(s in 1:S){

    ### step1: K by Metropolis-within-Gibbs

        # use random walk for K 7:3 Cauchy:uniform                                                                                                                                       
        ss <- runif(1)
        if (ss < 0.7) {  
          jump <- sample(-1:1, 1, prob = rep(1 / 3, 3))
        }
        else {
          jump <- round(rt(1, 1))  #discrete Cauchy
        }
        K.star <- K + jump
        while (K.star < (d + 2) || K.star > Kmax) {  # A bit hacky to ensure k doesn't go out of bounds
          if (ss < 0.7) {
            jump <- sample(-1:1, 1, prob = rep(1 / 3, 3))  # Ordinary proposal
          }
          else {
            jump <- round(rt(1, 1))  # Bold proposal
          }
          K.star <- K + jump
        }

        if(K.star!=K){# the weights will change
            G.star=rep(0,K.star)
            DPW.star=weights2(K.star,L,V,Z)
            G.star= DPW.star$G

            ### knots.star
            DPK.star=weights2(K.star-d,L,U,C)
            knot_diff.star=c(0,DPK.star$G)
            in_knots.star=rep(0,(K.star-d+1))
            in_knots.star=cumsum(knot_diff.star)

            ### the Bspline mixture density
            Bspline_dens=matrix(nr=K,nc=n)
            Bspline_dens.star=matrix(nr=K.star,nc=n)

            ### calculate the log posterior 
            XB=X%*%beta
            s1=s2=0
            bspline_dens=bspline_dens.star=NULL
            for(i in 1:n){
                w=exp(XB[i])/(1+exp(XB[i]))
                bspline_dens=dbspline(w,in_knots)
                Bspline_dens[,i]=bspline_dens
                bspline_dens.star=dbspline(w,in_knots.star)
                Bspline_dens.star[,i]=bspline_dens.star
                s1=s1+dnorm(Y[i],mu_f(t(G.star)%*%bspline_dens.star,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                s2=s2+dnorm(Y[i],mu_f(t(G)%*%bspline_dens,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
            }

            ### the accept or reject for K.star
            lnr=s1-s2+(K.star-K)*log(lambda)+log(factorial(K))-log(factorial(K.star))
            if(min(exp(lnr),1)>runif(1)){
                # store the outcome
                K=K.star
                G=NULL
                G=G.star
                DPW=DPK=NULL
                DPW=DPW.star
                DPK=DPK.star
                in_knots=NULL
                in_knots=in_knots.star
                Bspline_dens=matrix(nr=K.star,nc=n)
                Bspline_dens=Bspline_dens.star
                }
            # the adjustment of the dk
            if(exp(lnr)>0.65) dk=min(4,dk*2)
            else if(exp(lnr)<0.15) dk=max(2,dk/2)
        } # K.star!=K
        ### outcome store
        KK[s]=K 
    

    ### step2: lambda sample directly

        lambda=rgamma(1,a+K,b+1)
        Lambda[s]=lambda

  
    ### step3: r by Metropolis-within-Gibbs (q is sampled directly)

        # update rj in random order
        for(j in sample(1:p,p)){

            # sample q directly
            q=rbeta(1,1+r[j],2-r[j]) 
            rj.star=rbinom(1,1,q) # the proposal
            r.star=NULL
            r.star=r
            r.star[j]=rj.star

            ### if Bspline mixture density is NULL
            if(length(bspline_dens)==0){ 
                Bspline_dens=matrix(nr=K,nc=n)
                XB=X%*%beta
                for(i in 1:n){
                w=exp(XB[i])/(1+exp(XB[i]))
                bspline_dens=dbspline(w,in_knots)
                Bspline_dens[,i]=bspline_dens
                }
            }
    
            if(rj.star!=r[j] & sum(r.star)>=1){
                # updating q
                q.star=rbeta(1,1+rj.star,2-rj.star)
                # updating beta accordingly
                beta.star=NULL
                if(sum(r.star)==1){# one dimension
                    beta.star=rep(0,p)
                    beta.star[which(r.star==1)]=1
                    }
                else{
                    beta.star=beta
                    non_zeron=which(Beta[1:(s-1),j]!=0)
                    if(length(non_zeron)<=1){
                        beta.star[j]=rj.star*rnorm(1,0,1)
                    }
                    else{
                        beta.star[j]=rj.star*rnorm(1,mean(Beta[non_zeron,j]),sd(Beta[non_zeron,j]))
                    }
                # the normalisation
                beta.star=beta.star/c(sqrt(crossprod(beta.star)));
                if(beta.star[which(r.star==1)[1]]<0){beta.star=-beta.star}
                }# sum(r.star)>1
                nr=sum(r)
                nr.star=sum(r.star)

                ### calculate the log posterior 
                s1=s2=0
                Bspline_dens.star=matrix(nr=K,nc=n)
                XB.s=X%*%beta.star
                bspline_dens.star=NULL
                for(i in 1:n){
 		            w.star=exp(XB.s[i])/(1+exp(XB.s[i]))
                    bspline_dens.star=dbspline(w.star,in_knots)
                    Bspline_dens.star[,i]=bspline_dens.star

 			        s1=s1+dnorm(Y[i],mu_f(t(G)%*%bspline_dens.star,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 			        s2=s2+dnorm(Y[i],mu_f(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 		        }   

                ### the accept or reject for r.star
                lnr=s1-s2+(rj.star)*log(q.star/(1-q.star))-(r[j])*log(q/(1-q))+log(1-q.star)-log(1-q)+(nr-nr.star)*log(pi)/2+log(gamma(nr.star/2))-log(gamma(nr/2))
                if(min(exp(lnr),1)>runif(1)){
                    r[j]=rj.star
                    beta=NULL
                    beta=beta.star
                    Bspline_dens=matrix(nr=K,nc=n)
                    Bspline_dens=Bspline_dens.star
                }
            }# rj.star!=r[j] & sum(r.star)>=1
        } # j ends
        ### outcome store
        R[s,]=r
  
  
    ### step4: beta by Metropolis-within-Gibbs (the rho is updated; rho is the multivariate normal distribution of non zero beta)

        non_zerob=which(r==1)
        if(length(non_zerob)>1){

            # the beta updating
            beta.star=rep(0,p)
            # the rho updating
            mcmc_beta_r_rho=max(0,mcmc_beta_r_rho+(0.234-mcmc_r_beta_ratio)/sqrt(s))
            Rho[s]=mcmc_beta_r_rho

            # the multivariate normal distribution
            beta.star[non_zerob]=as.vector(rmvnorm(1,sqrt(2)*mcmc_beta_r_rho*beta[non_zerob],diag(rep(length(non_zerob)))))
            beta.star=beta.star/sqrt(c(crossprod(beta.star[non_zerob])))
            if(beta.star[which(beta.star!=0)[1]]<0){beta.star=-beta.star}

            Bspline_dens.star=matrix(nr=K,nc=n)
            XB.s=X%*%beta.star

            ### calculate the log posterior 
            s1=s2=0
            bspline_dens.star=NULL
 		    for(i in 1:n){
 		        w.star=exp(XB.s[i])/(1+exp(XB.s[i]))
                bspline_dens.star=dbspline(w.star,in_knots)
                Bspline_dens.star[,i]=bspline_dens.star
 			    s1=s1+dnorm(Y[i],mu_f(t(G)%*%bspline_dens.star,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 			    s2=s2+dnorm(Y[i],mu_f(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 		    }   
            
            ### the accept or reject for beta.star
            lnr=s1-s2+(t(beta.star-sqrt(2)*mcmc_beta_r_rho*beta)%*%diag(rep(p))%*%(beta.star-sqrt(2)*mcmc_beta_r_rho*beta)-t(beta-sqrt(2)*mcmc_beta_r_rho*beta.star)%*%diag(rep(p))%*%(beta-sqrt(2)*mcmc_beta_r_rho*beta.star))/2
            # updating ratio immediately
            mcmc_r_beta_ratio=min(exp(lnr),1)
            R_beta_ratio[s]=mcmc_r_beta_ratio

            if(mcmc_r_beta_ratio>runif(1)){
                beta=NULL;beta=beta.star;
                Bspline_dens=matrix(nr=K,nc=n)
                Bspline_dens=Bspline_dens.star
            }
        }# length(non_zerob)>1
        Beta[s,]=beta
  
    
    ### step5: V by Metropolis-within-Gibbs (V[l] l is randomly updated)

        for(l in sample(1:L,L)){

            # random walk of V[]
            dv=l/(l+2*sqrt(n))
            Lv=max(0,V[l]-dv)
            Uv=min(1,V[l]+dv)
            Vl.star=runif(1,Lv,Uv)
            V.star=V
            V.star[l]=Vl.star

            # G is updated
            G.star=NULL
            DPW.star=weights2(K,L,V.star,Z)
            G.star= DPW.star$G

            # the log posterior 
            s1=s2=0
            for(i in 1:n){
                s1=s1+dnorm(Y[i],mu_f(t(G.star)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                s2=s2+dnorm(Y[i],mu_f(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
            }

            ### the accept or reject for V.star
            lnr=s1-s2+(M-1)*(log(1-Vl.star)-log(1-V[l]))
            if(min(exp(lnr),1)>runif(1)){
                V[l]=Vl.star
                G=NULL
                G=G.star
                DPW=NULL
                DPW=DPW.star
            }
        }# l ends
        VV[s,]=V
  
    ### step6: Z by Metropolis-within-Gibbs (Z[l] l is randomly updated)

        for(l in sample(1:(L+1),L+1)){

            dz=l/(l+2*sqrt(n))
            Lz=max(0,Z[l]-dz)
            Uz=min(1,Z[l]+dz)
            Zl.star=runif(1,Lz,Uz)

            G.star=NULL
            G.star=G
            Zg=DPW$Zg
            Zg.star=Zg

            for (i in 1:K) { #if the G will change
                if ((Zl.star>(i-1)/K) & (Zl.star<=i/K)){
                Zg.star[l]=i
                break
                }
            }
            if(Zg.star[l]!=Zg[l]){

                # updating weights
                PW=DPW$P
                G.star[Zg[l]]=G.star[Zg[l]]-PW[l]
                G.star[Zg.star[l]]=G.star[Zg.star[l]]+PW[l]

                s1=s2=0
                for (i in 1:n){
                    s1=s1+dnorm(Y[i],mu_f(t(G.star)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                    s2=s2+dnorm(Y[i],mu_f(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                }
                lnr=s1-s2
                if(min(exp(lnr),1)>runif(1)){
                    Z[l]=Zl.star
                    G=G.star
                    DPW=NULL
                    DPW=list(G=G.star,Zg=Zg.star,P=PW)}
            }
        }# l ends
        ZZ[s,]=Z
        # with V,Z ready G.star is stored
        GG[s,1:K]=G


    ### step7: U by Metropolis-within-Gibbs (The knots updates))

        for(l in sample(1:L,L)){

 	        du=l/(l+2*sqrt(n))
            Lu=max(0,U[l]-du)
            Uu=min(1,U[l]+du)
 	        Ul.star=runif(1,Lu,Uu)
            U.star=U;U.star[l]=Ul.star;

            ### knots.star
            DPK.star=weights2(K-d,L,U.star,C)
            knot_diff.star=c(0,DPK.star$G)
            in_knots.star=rep(0,(K-d+1))
            in_knots.star=cumsum(knot_diff.star)

            Bspline_dens.star=matrix(nr=K,nc=n)
            XB=X%*%beta

            s1=s2=0
            bspline_dens=bspline_dens.star=NULL
            for(i in 1:n){
                w=exp(XB[i])/(1+exp(XB[i]))
                bspline_dens.star=dbspline(w,in_knots.star)
                Bspline_dens.star[,i]=bspline_dens.star
                s1=s1+dnorm(Y[i],mu_f(t(G)%*%bspline_dens.star,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                s2=s2+dnorm(Y[i],mu_f(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
            }
 	        lnr=s1-s2+(M-1)*(log(1-Ul.star)-log(1-U[l]))                                                                                                                                                  
 	        if(min(exp(lnr),1)>runif(1)){
                U[l]=Ul.star
                in_knots=NULL
                in_knots=in_knots.star
                DPK=NULL
                DPK=DPK.star
                Bspline_dens=matrix(nr=K,nc=n)
                Bspline_dens=Bspline_dens.star}
        }## l ends
        UU[s,]=U


    ### step8: C by Metropolis-within-Gibbs (The knots updates))

        for(l in sample(1:(L+1),L+1)){

            dc=l/(l+2*sqrt(n))
            Lc=max(0,C[l]-dc)
            Uc=min(1,C[l]+dc)
            Cl.star=runif(1,Lc,Uc)

            #knots.star
            in_knots.star=NULL
            in_knots.star=in_knots
            Cg=DPK$Zg
            Cg.star=Cg
            for (i in 1:(K-d)) { # if the knots will change
                if ((Cl.star>(i-1)/(K-d)) & (Cl.star<=i/(K-d))){
                Cg.star[l]=i
                break
                }
            }
            if(Cg.star[l]!=Cg[l]){
                PK=DPK$P
                GK.s=DPK$G
                GK.s[Cg[l]]=GK.s[Cg[l]]-PK[l]
                GK.s[Cg.star[l]]=GK.s[Cg.star[l]]+PK[l]

                knot_diff.star=c(0,GK.s)
                in_knots.star=rep(0,(K-d+1))
                in_knots.star=cumsum(knot_diff.star)
                Bspline_dens.star=matrix(nr=K,nc=n)
                XB=X%*%beta
                s1=s2=0
                bspline_dens=bspline_dens.star=NULL

                for(i in 1:n){
                    w=exp(XB[i])/(1+exp(XB[i]))
                    bspline_dens.star=dbspline(w,in_knots.star)
                    Bspline_dens.star[,i]=bspline_dens.star
                    s1=s1+dnorm(Y[i],mu_f(t(G)%*%bspline_dens.star,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                    s2=s2+dnorm(Y[i],mu_f(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
                }
 		        lnr=s1-s2
            if(min(exp(lnr),1)>runif(1)){
                C[l]=Cl.star
                in_knots=NULL
                in_knots=in_knots.star
                DPK=NULL
                DPK=list(G=GK.s,Zg=Cg.star,P=PK)
                Bspline_dens=matrix(nr=K,nc=n)
                Bspline_dens=Bspline_dens.star}
            }
        } # l ends
        CC[s,]=C
        GG_knots[s,1:(K-d)]=t(DPK$G)

    
    ### step9: sigma& kc (sampled directly)
    
        mu=t(Y)-mu_f(t(G)%*%Bspline_dens,baseb)
        sye=mu%*%diag(1/kc)%*%t(mu)
        sigma=rinvgamma(1,(nu0+3*n)/2,(nu0*s0+2*sum(kc)+sye/8)/2)
        Sigma[s]=sigma
  

        eta2=2/sigma
        for (i in 1:n){
        eta1=((Y[i]-mu_f(t(G)%*%Bspline_dens[,i],baseb))^2)/(8*sigma)
            kc[i]=rgig(1,1/2,sqrt(eta1),sqrt(eta2))
        } 
        KC[s,]=kc
  
     print(paste(s,Sys.time()," "))
    # END: MCMC loop
    }

    ###
    # Which iterations to keep
    ###

    keep=seq(burnin + 1, S)
    KK=KK[keep]
    R=R[keep,]
    Beta=Beta[keep,]
    VV=VV[keep,]
    ZZ=ZZ[keep,]
    UU=UU[keep,]
    CC=CC[keep,]
    GG=GG[keep,]
    GG_knots=GG_knots[keep,]
    Lambda=Lambda[keep]
    Sigma=Sigma[keep]
    KC=KC[keep,]

    
    ###
    # the estimate
    ###

    R.pred=Beta.pred=G.est=Gknots.est=NULL
    ### the function to get mode
    getmode = function(v) {
        uniqv = unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    # the median of KK is the point estimate
    K.pred=getmode(KK);

    # the R.pred is also the median
    nzk=which(KK==K.pred)
    lnzk=length(nzk)
    R0=R[nzk,]
    Beta0=Beta[nzk,]
    Sigma0=Sigma[nzk]
    R.uniqv=unique(R0)
    lR=nrow(R.uniqv)
    R.count=rep(0,lR)
    R.cg=rep(0,lnzk)
    for (i in 1: lnzk){
        for (j in 1: lR){
            if (isTRUE(all.equal(R0[i,],R.uniqv[j,]))){
                R.count[j]=R.count[j]+1
                R.cg[i]=j
            }
        }
    } 
    R.pred=R.uniqv[which.max(R.count),]

    # the beta.pred
    beta0.order=which(R.cg==which.max(R.count))
    Beta00=matrix(nr=max(R.count),nc=p)
    Beta00=Beta0[beta0.order,]
    beta=apply(Beta00,2,mean,na.rm=TRUE)
    beta=beta/c(sqrt(crossprod(beta)))
    if(beta[which(R.pred==1)[1]]<0){beta=-beta}
    Beta.est=beta

    # the G
    G.est[1:K.pred]=apply(GG[nzk,1:K.pred],2,mean,na.rm=TRUE)
    
    # the knots
    if(K.pred==4){
        in_knots=c(0,1)
        }
        else{
            Gknots.est[1:(K.pred-d)]=apply(GG_knots[nzk,1:(K.pred-d)],2,mean,na.rm=TRUE)
            }
    knot_diff.est=c(0,Gknots.est)
    In_knots.est=rep(0,(K.pred-d+1))
    In_knots.est=cumsum(knot_diff.est)

    # the KC sigma 
    kc=apply(KC[nzk[floor(lnzk/3):lnzk],],2,mean)
    Sigma00=rep(0,max(R.count))
    Sigma00=Sigma0[beta0.order]
    Sigma.est=mean(Sigma00)


    ###
    # the output
    ###

    out=list(
        K.est=K.pred,
        R.est=R.pred,
        Beta.est=Beta.est,
        G.est=G.est,
        In_knots.est=In_knots.est,
        Sigma.est=Sigma.est)

    return(out)

}







