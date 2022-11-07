# All competitors considered in the paper. TVD-OG-MBD-ED
library(fda)
library(roahd)
library(robustbase)
library(fdaoutlier)
library(mrfDepth)
set.seed(167)
############################################################################################
TE=500 #Total experiments
met=4
TPR2=matrix(NaN,ncol = TE,nrow = met )
FPR2=matrix(NaN,ncol = TE,nrow = met)
TPR22=matrix(NaN,ncol = TE,nrow = met )
FPR22=matrix(NaN,ncol = TE,nrow = met)
TPR23=matrix(NaN,ncol = TE,nrow = met )
FPR23=matrix(NaN,ncol = TE,nrow = met)
TPR24=matrix(NaN,ncol = TE,nrow = met )
FPR24=matrix(NaN,ncol = TE,nrow = met)
TPR25=matrix(NaN,ncol = TE,nrow = met )
FPR25=matrix(NaN,ncol = TE,nrow = met)
TPR26=matrix(NaN,ncol = TE,nrow = met )
FPR26=matrix(NaN,ncol = TE,nrow = met)
TPR27=matrix(NaN,ncol = TE,nrow = met )
FPR27=matrix(NaN,ncol = TE,nrow = met)
M=100 #number of considered time points
N=1000 # Number of curves
k=6
#Generations of the models  
cov.fun=function(d,k,c,mu){
  k*exp(-(1/c)*d^mu)
}
t=seq(0,1,len=M)
d1=as.matrix(dist(t,upper=TRUE,diag=TRUE))
c=1
sig=1
for (u in 1:TE) {
  #Generate the models 
  D=rbinom(N,1,0.1)
  # parameters of the covariance function 
  c=1
  sig=1
  t.cov=cov.fun(d1,sig,c,1)
  # Cholesky Decomposition
  L=chol(t.cov)
  e=matrix(rnorm(N*M),M,N)
  #Model_0 Base Model
  xt1=as.matrix(t(L)%*%e)
  # Model 1. Pure Dependnece Shape Outliers
  #covariance structure of dependence shape outliers
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
  #model 2 Phase shape outliers
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
  
  #model 3 High frequency low amplitude
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
  #Model 4  slightly modification of the Arribas-Romo model
  xt9=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      xt9[,i] =30*t*(1-t)^(3/2)+t(L2)%*%e[,i]
    }else{
      xt9[,i] =30*t*(1-t)^(3/2)+t(L)%*%e[,i]
    }
  }

  #model 5 Central high frequency-low amplitude
  xt10=matrix(NaN,M,N)
  for (i in 1:N) {
    if(D[i]==1){
      theta=runif(1,0.25,0.5)
      xt10[,i] =0.7*sin(40*(t+theta)*pi)+t(L2)%*%e[,i]
    }else{
      xt10[,i] =xt1[,i]
    }
  }

  model=list(xt5,xt7,xt8,xt9,xt10)
  E=length(model)
  #other methods
  other_methods_outliers=list()
  # #1. TVD
  other_methods_outliers[[1]]=tvdmss(t(model[[1]]))$outliers
  # #2. MBD
  other_methods_outliers[[2]]=fda::fbplot(model[[1]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers[[3]]=functional_boxplot(t(model[[1]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[1]]))
  other_methods_outliers[[4]]=outliergram(fD,display = F,)$ID_outliers

  EE=length(other_methods_outliers)
  true_out2=rep(0,EE)
  false_out2=rep(0,EE)
  C=D
  h = 1
  for (i in other_methods_outliers) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out2[h] = true_out2[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out2[h] =false_out2[h] + 1
      }
    }
    h = h + 1
  }
  
  total_out=sum(C)
  TPR2[,u]=(true_out2/total_out)*100
  FPR2[,u]=(false_out2/(N-total_out))*100 
  #########################################
  
  other_methods_outliers2=list()
  
  # #1. TVD+MSS
  other_methods_outliers2[[1]]=tvdmss(t(model[[2]]))$outliers
  # #2. MBD
  other_methods_outliers2[[2]]=fda::fbplot(model[[2]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers2[[3]]=functional_boxplot(t(model[[2]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[2]]))
  other_methods_outliers2[[4]]=outliergram(fD,display = F,)$ID_outliers

  true_out22=rep(0,EE)
  false_out22=rep(0,EE)
  
  h = 1
  for (i in other_methods_outliers2) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out22[h] = true_out22[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out22[h] =false_out22[h] + 1
      }
    }
    h = h + 1
  }
  
  total_out=sum(C)
  TPR22[,u]=(true_out22/total_out)*100
  FPR22[,u]=(false_out22/(N-total_out))*100 
  ##############################################
  other_methods_outliers3=list()
  
  # #1. TVD+MSS
  other_methods_outliers3[[1]]=tvdmss(t(model[[3]]))$outliers
  # #2. MBD
  other_methods_outliers3[[2]]=fda::fbplot(model[[3]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers3[[3]]=functional_boxplot(t(model[[3]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[3]]))
  other_methods_outliers3[[4]]=outliergram(fD,display = F,)$ID_outliers
  true_out23=rep(0,EE)
  false_out23=rep(0,EE)
  
  h = 1
  for (i in other_methods_outliers3) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out23[h] = true_out23[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out23[h] =false_out23[h] + 1
      }
    }
    h = h + 1
  }
  
  total_out=sum(C)
  TPR23[,u]=(true_out23/total_out)*100
  FPR23[,u]=(false_out23/(N-total_out))*100 
  #############################################################
  other_methods_outliers4=list()
  
  # #1. TVD+MSS
  other_methods_outliers4[[1]]=tvdmss(t(model[[4]]))$outliers
  # #2. MBD
  other_methods_outliers4[[2]]=fda::fbplot(model[[4]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers4[[3]]=functional_boxplot(t(model[[4]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[4]]))
  other_methods_outliers4[[4]]=outliergram(fD,display = F,)$ID_outliers

  true_out24=rep(0,EE)
  false_out24=rep(0,EE)
  
  h = 1
  for (i in other_methods_outliers4) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out24[h] = true_out24[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out24[h] =false_out24[h] + 1
      }
    }
    h = h + 1
  }
  
  total_out=sum(C)
  TPR24[,u]=(true_out24/total_out)*100
  FPR24[,u]=(false_out24/(N-total_out))*100 
  
  #############################################################
  other_methods_outliers5=list()
  
  # #1. TVD+MSS
  other_methods_outliers5[[1]]=tvdmss(t(model[[5]]))$outliers
  # #2. MBD
  other_methods_outliers5[[2]]=fda::fbplot(model[[5]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers5[[3]]=functional_boxplot(t(model[[5]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[5]]))
  other_methods_outliers5[[4]]=outliergram(fD,display = F,)$ID_outliers

  true_out25=rep(0,EE)
  false_out25=rep(0,EE)
  
  h = 1
  for (i in other_methods_outliers5) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out25[h] = true_out25[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out25[h] =false_out25[h] + 1
      }
    }
    h = h + 1
  }
  total_out=sum(C)
  TPR25[,u]=(true_out25/total_out)*100
  FPR25[,u]=(false_out25/(N-total_out))*100 
  
  #############################################################
  other_methods_outliers6=list()
  
  # #1. TVD+MSS
  other_methods_outliers6[[1]]=tvdmss(t(model[[6]]))$outliers
  # #2. MBD
  other_methods_outliers6[[2]]=fda::fbplot(model[[6]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers6[[3]]=functional_boxplot(t(model[[6]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[6]]))
  other_methods_outliers6[[4]]=outliergram(fD,display = F,)$ID_outliers

  true_out26=rep(0,EE)
  false_out26=rep(0,EE)
  
  h = 1
  for (i in other_methods_outliers6) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out26[h] = true_out26[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out26[h] =false_out26[h] + 1
      }
    }
    h = h + 1
  }
  
  total_out=sum(C)
  TPR26[,u]=(true_out26/total_out)*100
  FPR26[,u]=(false_out26/(N-total_out))*100 
  
  #############################################################
  other_methods_outliers7=list()
  
  # #1. TVD+MSS
  other_methods_outliers7[[1]]=tvdmss(t(model[[7]]))$outliers
  # #2. MBD
  other_methods_outliers7[[2]]=fda::fbplot(model[[7]],plot = F,method = 'MBD')$outpoint
  # #3. ED
  other_methods_outliers7[[3]]=functional_boxplot(t(model[[7]]),depth_method = 'extremal',central_region = 0.5)$outliers
  # #4. OG+MBD
  fD=fData(t,t(model[[7]]))
  other_methods_outliers7[[4]]=outliergram(fD,display = F,)$ID_outliers

  true_out27=rep(0,EE)
  false_out27=rep(0,EE)
  
  h = 1
  for (i in other_methods_outliers7) {
    for (j in i) {
      if (j %in% which(C == 1)) {
        true_out27[h] = true_out27[h] + 1
      } else if (j %in% which(C == 0)) {
        false_out27[h] =false_out27[h] + 1
      }
    }
    h = h + 1
  }
  total_out=sum(C)
  TPR27[,u]=(true_out27/total_out)*100
  FPR27[,u]=(false_out27/(N-total_out))*100 

}

#Other models
#TPR
mean.TPR2=apply(TPR2,1,mean);mean.TPR2
sd.TPR2=apply(TPR2,1,sd);sd.TPR2
#False positive rate
mean.FPR2=apply(FPR2,1,mean);mean.FPR2
sd.FPR2=apply(FPR2,1,sd);sd.FPR2


#TPR
mean.TPR22=apply(TPR22,1,mean);mean.TPR22
sd.TPR22=apply(TPR22,1,sd);sd.TPR22
#False positive rate
mean.FPR22=apply(FPR22,1,mean);mean.FPR22
sd.FPR22=apply(FPR22,1,sd);sd.FPR22

#TPR
mean.TPR23=apply(TPR23,1,mean);mean.TPR23
sd.TPR23=apply(TPR23,1,sd);sd.TPR23
#False positive rate
mean.FPR23=apply(FPR23,1,mean);mean.FPR23
sd.FPR23=apply(FPR23,1,sd);sd.FPR23

#TPR
mean.TPR24=apply(TPR24,1,mean);mean.TPR24
sd.TPR24=apply(TPR24,1,sd);sd.TPR24
#False positive rate
mean.FPR24=apply(FPR24,1,mean);mean.FPR24
sd.FPR24=apply(FPR24,1,sd);sd.FPR24


#TPR
mean.TPR25=apply(TPR25,1,mean);mean.TPR25
sd.TPR25=apply(TPR25,1,sd);sd.TPR25
#False positive rate
mean.FPR25=apply(FPR25,1,mean);mean.FPR25
sd.FPR25=apply(FPR25,1,sd);sd.FPR25

#TPR
mean.TPR26=apply(TPR26,1,mean);mean.TPR26
sd.TPR26=apply(TPR26,1,sd);sd.TPR26
#False positive rate
mean.FPR26=apply(FPR26,1,mean);mean.FPR26
sd.FPR26=apply(FPR26,1,sd);sd.FPR26

#TPR
mean.TPR27=apply(TPR27,1,mean);mean.TPR27
sd.TPR27=apply(TPR27,1,sd);sd.TPR27
#False positive rate
mean.FPR27=apply(FPR27,1,mean);mean.FPR27
sd.FPR27=apply(FPR27,1,sd);sd.FPR27

