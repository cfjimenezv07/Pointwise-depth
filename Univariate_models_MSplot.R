#Univariate models compared with MS plot

options(scipen=100)
set.seed(167)
TE=500
M=100
N=1000
k=6
cov.fun=function(d,k,c,mu){
  k*exp(-(1/c)*d^mu)
}
combinat=function(n,p){if (n<p){combinat=0}else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}}
t=seq(0,1,len=M)
d1=as.matrix(dist(t,upper=TRUE,diag=TRUE))
TPR=matrix(NaN,ncol = TE,nrow =5)
FPR=matrix(NaN,ncol = TE,nrow =5)
for (u in 1:TE) {
  #Generations of the models 
  D=rbinom(N,1,0.1)
  #covariance function in time
  c=1
  sig=1
  t.cov=cov.fun(d1,sig,c,1)
  # Cholesky Decomposition
  L=chol(t.cov)
  mu=0
  e=matrix(rnorm(N*M),M,N)
  #Model_0 Base Model
  xt1=as.matrix(mu+t(L)%*%e)
  #matplot(xt1,type='l',col=9,ylim = c(-8,8))
  
  #Model 1. Pure Shape Outliers
  #covariance structure of outliers
  t.cov1=cov.fun(d1,k,c,0.1)
  # Cholesky Decomposition
  L1=chol(t.cov1)
  out5=as.matrix(t(L1)%*%e)
  xt5=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      xt5[,i]=out5[,i]
    }else{
      xt5[,i]=xt1[,i]
    }
  }
  mm=matrix(t,M,N)
  
  
  
  
  
  #model 2
  m_base=function(t){2*sin(15*pi*t)}
  m_cont=function(t){2*sin(15*pi*t+4)}
  out7.0=apply(mm,2,m_base)+t(L)%*%e
  out7.1=apply(mm,2,m_cont)+t(L)%*%e
  xt7=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      xt7[,i] =out7.1[,i]
    }else{
      xt7[,i] =out7.0[,i]
    }
  }
  
  
  
  
  #model 3
  A=0.1
  B=1
  m_g=function(t){A+B*atan(t)}
  m_o=function(t){atan(t)}
  t.cov2=cov.fun(d1,0.1,4,0.1)
  # Cholesky Decomposition
  L2=chol(t.cov2)
  out8.0=apply(mm,2,m_g)+t(L)%*%e
  out8.1=apply(mm,2,m_o)+t(L2)%*%e
  xt8=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      xt8[,i] =out8.1[,i]
    }else{
      xt8[,i] =out8.0[,i]
    }
    
  }
  
  
  
  
  
  #Model 4
  xt9=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      xt9[,i] =30*t*(1-t)^(3/2)+t(L2)%*%e[,i]
    }else{
      xt9[,i] =30*t*(1-t)^(3/2)+t(L)%*%e[,i]
    }
  }
  
  
  
  
  #model 5
  xt10=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      theta=runif(1,0.25,0.5)
      xt10[,i] =0.1*sin(40*(t+theta)*pi)+t(L2)%*%e[,i]
    }else{
      xt10[,i] =xt1[,i]
    }
  }
  
  
  model=list(xt5,xt7,xt8,xt9,xt10)
  E=length(model)
  
  
  
  outliers=list()
  for (y in 1:E) {
    msp= fdaoutlier::msplot(t(model[[y]]),plot=F)
    outliers[[y]]= msp$outliers
  }
  
  true_out=rep(0,E)
  false_out=rep(0,E)
  
  
  
  o = 1
  for (i in outliers) {
    for (j in i) {
      if (j %in% which(D == 1)) {
        true_out[o] = true_out[o] + 1
      } else if (j %in% which(D == 0)) {
        false_out[o] =false_out[o] + 1
      }
    }
    o = o + 1
  }
  
  
  total_out=sum(D)
  
  TPR[,u]=(true_out/total_out)*100
  FPR[,u]=(false_out/(N-total_out))*100
  
}


# True positive rate
mean.TPR=apply(TPR,1,mean);mean.TPR
sd.TPR=apply(TPR,1,sd);sd.TPR

# False positive rate
mean.FPR=apply(FPR,1,mean);mean.FPR
sd.FPR=apply(FPR,1,sd);sd.FPR

