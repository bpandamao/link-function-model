
library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma)

##### data structure
gibbs_bernstein_est_linkfunction= function(
  X,Y,S=10000,L=max(20,n^(1/3)),Kmax=50,M=1,a=0.01,b=0.01,K=10,nu0=500,s0=0.05,lambda=10,dk=4,sigma=0.5,baseb=1.15,burnin=2000,logstate=TRUE
){
n=nrow(X)
p=ncol(X)

##### The initialization

# Initialise G
G=rep(NA,K);
# Initialise V,Z
V = rep(1/2,L)
Z = runif(L + 1)

#obtain G(j)
DP=weights2(K,L,V,Z)
G=DP$G

#Initialise r
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
    return(mu_f(x,baseb))
  } 
}else{
  mu_f=function(x,baseb){
    return(x)
  }
}

##### Define Parameter space
KK=Lambda=beta_dens=Raa=Sigma=AA=NULL
R=Beta=matrix(nr=S,nc=p)
VV=matrix(nr=S,nc=L)
ZZ=matrix(nr=S,nc=L+1)
KC=matrix(nr=S,nc=n)
GG=matrix(nr=S,nc=Kmax)
Rho=R_beta_ratio=NULL

##### the updating process

for(s in 1:S){
### the K
 Lk=max(3,round(K-dk));Uk=min(Kmax,ceiling(K+dk));K.s=sample(Lk:Uk,1);
 if(K.s!=K){
   G.s=rep(0,K.s)
   DP.s=weights2(K.s,L,V,Z)
   G.s= DP.s$G

 Beta_dens.s=matrix(nr=K.s,nc=n);Beta_dens=matrix(nr=K,nc=n)
 XB=X%*%beta
 s1=s2=0
 beta_dens=beta_dens.s=NULL
 for(i in 1:n){
   w=exp(XB[i])/(1+exp(XB[i]))
   for (j in 1:K.s){
     beta_dens.s[j]=dbeta(w,j,K.s-j+1)
   }
    Beta_dens.s[,i]=beta_dens.s
   for (g in 1:K){
     beta_dens[g]=dbeta(w,g,K-g+1)
   }
    Beta_dens[,i]=beta_dens
    s1=s1+dnorm(Y[i],mu_f(t(G.s)%*%beta_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
    s2=s2+dnorm(Y[i],mu_f(t(G)%*%beta_dens,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 }
 
 lnr=s1-s2+(K.s-K)*log(lambda)+log(factorial(K))-log(factorial(K.s))
 if(min(exp(lnr),1)>runif(1)){
   K=K.s
   G=NULL
   G=G.s
   DP=NULL
   DP=DP.s
   Beta_dens=matrix(nr=K.s,nc=n)
   Beta_dens=Beta_dens.s}
 if(exp(lnr)>0.65) dk=min(4,dk*2) 
 else if(exp(lnr)<0.15) dk=max(2,dk/2)
 } ### K.s!=K
 KK[s]=K

### the lambda
lambda=rgamma(1,a+K,b+1)
Lambda[s]=lambda


 ### the r
 for(j in sample(1:p,p)){
 	q=rbeta(1,1+r[j],2-r[j]) 
   rj.s=rbinom(1,1,q) #0,1
   r.s=NULL;r.s=r
   r.s[j]=rj.s
    #### if Beta_dens is NULL
   if(length(beta_dens)==0){ 
        Beta_dens=matrix(nr=K,nc=n);
        XB=X%*%beta
        for(i in 1:n){
          w=exp(XB[i])/(1+exp(XB[i]))
          for (jj in 1:K){
          beta_dens[jj]=dbeta(w,jj,K-jj+1)
          }
          Beta_dens[,i]=beta_dens
        }
   }
 	if(rj.s!=r[j] & sum(r.s)>=1){
 	  q.s=rbeta(1,1+rj.s,2-rj.s)
 		beta.s=NULL
 		if(sum(r.s)==1){
       beta.s=rep(0,p)
       beta.s[which(r.s==1)]=1
       }
       else{
         beta.s=beta
         wns=which(Beta[1:(s-1),j]!=0)
         if(length(wns)<=1){
           beta.s[j]=rj.s*rnorm(1,0,1)
           }
           else{
             beta.s[j]=rj.s*rnorm(1,mean(Beta[wns,j]),sd(Beta[wns,j]))
             }
 			  beta.s=beta.s/c(sqrt(crossprod(beta.s)));
 			  if(beta.s[which(r.s==1)[1]]<0){beta.s=-beta.s};
        }### sum(r.s)>1
 		nr=sum(r);nr.s=sum(r.s)
    s1=s2=0;Beta_dens.s=matrix(nr=K,nc=n)
 		XB.s=X%*%beta.s
    beta_dens.s=NULL
 		for(i in 1:n){
 		  w.s=exp(XB.s[i])/(1+exp(XB.s[i]))
      for (jj in 1:K){
        beta_dens.s[jj]=dbeta(w.s,jj,K-jj+1)
        }
 			s1=s1+dnorm(Y[i],mu_f(t(G)%*%beta_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 			s2=s2+dnorm(Y[i],mu_f(t(G)%*%Beta_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 			Beta_dens.s[,i]=beta_dens.s
 		}   
 		lnr=s1-s2+(rj.s)*log(q.s/(1-q.s))-(r[j])*log(q/(1-q))+log(1-q.s)-log(1-q)+(nr-nr.s)*log(pi)/2+log(gamma(nr.s/2))-log(gamma(nr/2))
 		if(min(exp(lnr),1)>runif(1)){
       r[j]=rj.s
       beta=NULL;beta=beta.s
       Beta_dens=matrix(nr=K,nc=n);Beta_dens=Beta_dens.s
       }
 	}### rj.s!=r[j] & sum(r.s)>=1
 } ### j ends
 R[s,]=r


### the beta
 wb=which(r==1)
 if(length(wb)>1){
  beta.s=rep(0,p)
  mcmc_beta_r_rho=max(0,mcmc_beta_r_rho+(0.234-mcmc_r_beta_ratio)/sqrt(s))
  Rho[s]=mcmc_beta_r_rho
  beta.s[wb]=as.vector(rmvnorm(1,sqrt(2)*mcmc_beta_r_rho*beta[wb],diag(rep(length(wb)))))
  beta.s=beta.s/sqrt(c(crossprod(beta.s[wb])))
  if(beta.s[which(beta.s!=0)[1]]<0){beta.s=-beta.s}
  Beta_dens.s=matrix(nr=K,nc=n)
  XB.s=X%*%beta.s
  s1=s2=0
  beta_dens.s=NULL
 	for(i in 1:n){
 		w.s=exp(XB.s[i])/(1+exp(XB.s[i]))
    for (j in 1:K){
      beta_dens.s[j]=dbeta(w.s,j,K-j+1)
      }
 	  Beta_dens.s[,i]=beta_dens.s
 	  s1=s1+dnorm(Y[i],mu_f(t(G)%*%beta_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 	  s2=s2+dnorm(Y[i],mu_f(t(G)%*%Beta_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
    }
 lnr=s1-s2+(t(beta.s-sqrt(2)*mcmc_beta_r_rho*beta)%*%solve(diag(rep(p)))%*%(beta.s-sqrt(2)*mcmc_beta_r_rho*beta)-t(beta-sqrt(2)*mcmc_beta_r_rho*beta.s)%*%solve(diag(rep(p)))%*%(beta-sqrt(2)*mcmc_beta_r_rho*beta.s))/2
 mcmc_r_beta_ratio=min(exp(lnr),1);R_beta_ratio[s]=mcmc_r_beta_ratio
 if(mcmc_r_beta_ratio>runif(1)){
   beta=NULL;beta=beta.s;
   Beta_dens=matrix(nr=K,nc=n);Beta_dens=Beta_dens.s
   }
}##### length(wb)>1
 Beta[s,]=beta


### the V
 for(l in sample(1:L,L)){
 	dv=l/(l+2*sqrt(n));Lv=max(0,V[l]-dv);Uv=min(1,V[l]+dv)
 	Vl.s=runif(1,Lv,Uv);V.s=V;V.s[l]=Vl.s;G.s=NULL;
 	DP.s=weights2(K,L,V.s,Z)
 	G.s= DP.s$G
 	s1=s2=0
  for(i in 1:n){
 		s1=s1+dnorm(Y[i],mu_f(t(G.s)%*%Beta_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 		s2=s2+dnorm(Y[i],mu_f(t(G)%*%Beta_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 	}
 	lnr=s1-s2+(M-1)*(log(1-Vl.s)-log(1-V[l]))
 	if(min(exp(lnr),1)>runif(1)){V[l]=Vl.s;G=NULL;G=G.s;DP=NULL;DP=DP.s}
 }## l ends
 VV[s,]=V
  
 ### the Z
 for(l in sample(1:(L+1),L+1)){
 	dz=l/(l+2*sqrt(n));Lz=max(0,Z[l]-dz);Uz=min(1,Z[l]+dz)
 	Zl.s=runif(1,Lz,Uz);G.s=NULL;G.s=G;Zg=DP$Zg;Zg.s=Zg
  for (i in 1:K) {
    if ((Zl.s>(i-1)/K) & (Zl.s<=i/K)){
      Zg.s[l]=i
      break
    }
  }
  if(Zg.s[l]!=Zg[l]){
    P=DP$P
    G.s[Zg[l]]=G.s[Zg[l]]-P[l]
    G.s[Zg.s[l]]=G.s[Zg.s[l]]+P[l]
    s1=s2=0
    for (i in 1:n){
 		s1=s1+dnorm(Y[i],mu_f(t(G.s)%*%Beta_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 		s2=s2+dnorm(Y[i],mu_f(t(G)%*%Beta_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 	}
 		lnr=s1-s2
 		if(min(exp(lnr),1)>runif(1)){Z[l]=Zl.s;G=G.s;DP=NULL;DP=list(G=G,Zg=Zg.s,P=P)}
 }
}
 ZZ[s,]=Z;GG[s,1:K]=G

 ### the sigma
 mu=t(Y)-mu_f(t(G)%*%Beta_dens,baseb)
 sy=mu%*%diag(1/kc)%*%t(mu)
 sigma=rinvgamma(1,(nu0+3*n)/2,(nu0*s0+2*sum(kc)+sy/8)/2)
 Sigma[s]=sigma
   
 ### the kc
 eta2=2/sigma
  for (i in 1:n){
    eta1=((Y[i]-mu_f(t(G)%*%Beta_dens[,i],baseb))^2)/(8*sigma)
  kc[i]=rgig(1,1/2,sqrt(eta1),sqrt(eta2))
 } 
 KC[s,]=kc
 
 print(paste(s,Sys.time()," "))
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
GG=GG[keep,]
Lambda=Lambda[keep]
Sigma=Sigma[keep]
KC=KC[keep,]

###the posterior estimate
R.pred=Beta.pred=G.est=NULL

### the K.pred
### the function to get mode
getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
  }
K.pred=getmode(KK);
### the R.pred
wke=which(KK==K.pred)
lwke=length(wke)
R0=R[wke,]
Beta0=Beta[wke,]
R.uniqv=unique(R0)
lR=nrow(R.uniqv)
R.count=rep(0,lR)
R.cg=rep(0,lwke)
for (i in 1: lwke){
  for (j in 1: lR){
    if (isTRUE(all.equal(R0[i,],R.uniqv[j,]))){
      R.count[j]=R.count[j]+1
      R.cg[i]=j
    }
  }

} 
R.pred=R.uniqv[which.max(R.count),]

### the beta.pred
beta0.order=which(R.cg==which.max(R.count))
Beta00=matrix(nr=max(R.count),nc=p)
Beta00=Beta0[beta0.order,]
beta=apply(Beta00,2,mean,na.rm=TRUE)
beta=beta/c(sqrt(crossprod(beta)))
if(beta[which(R.pred==1)[1]]<0){beta=-beta}
Beta.est=beta

if(K.pred==1){G.est=1}else{G.est=apply(GG[wke,1:K.pred],2,mean,na.rm=TRUE)}

kc=apply(KC[wke[floor(lwke/3):lwke],],2,mean)
sigma.est=mean(Sigma[wke[floor(lwke/3):lwke]]);
Lambda.est=mean(Lambda)

out=list(
  K.est=K.pred,
  R.est=R.pred,
  Beta.est=Beta.est,
  G.est=G.est,
  Sigma.est=sigma.est)

return(out)


}
