#Multivariate models with MS plot
#####################################################################################################
library(fdaoutlier)
library(RandomFields)
library(ggplot2)
library(reshape2)
library(grid)
library(scatterplot3d)
library(pvclust)
library(MASS)
library(fda.usc)
library(ddalpha)
`%!in%`=Negate(`%in%`)
################################################################################
base_model = function (n, p, d, c){
  data = array(0, dim = c(n, p, d))
  z = seq(1, p) / p
  model = RMbiwm(nudiag = c(1.2, 0.6), nured = 1, rhored = 0.6, 
                 cdiag = c(1, 1) * 1, s = 10 * c(0.02, 0.016, 0.01))
  for (i in 1:n) {
    data[i, , ] = matrix(data = RFsimulate(model, z), nrow = p, ncol = 2)/10
  }
  return(data)
}
################################################################################
set.seed(167)
TE=500
TPR=matrix(NaN,ncol = TE,nrow =4)
FPR=matrix(NaN,ncol = TE,nrow = 4)
n = 100
p = 100
d = 2
c = 0.1
for (u in 1:TE) {
  
  #Base model
  data = base_model(n, p, d, c)
  z = seq(1, p) / p
  D=rbinom(n, 1, 0.1)
  o=which(D==1)
  A=data[,,1]
  B=data[,,2]
  #Contaminated models
  #1
  curves11=matrix(0,nrow = n,ncol = p)
  curves12=matrix(0,nrow = n,ncol = p)
  for (i in 1:n) {
    curves11[i,] =0.5*cos(80 * pi * z)/10+A[i,]
    curves12[i,] =0.75*sin(40 * pi * z)/10+B[i,]
  }
  c_model_1 = array(0,dim = c(n,p,2))
  c_model_1[,,1]=curves11
  c_model_1[,,2]=curves12
  
  #2
  curves21=matrix(0,nrow = n,ncol = p)
  curves22=matrix(0,nrow = n,ncol = p)
  for (i in 1:n) {
    curves21[i,] =2*cos(80 * pi * z)/10+A[i,]
    curves22[i,] =3*sin(40 * pi * z)/10+B[i,]
  }
  c_model_2 = array(0,dim = c(n,p,2))
  c_model_2[,,1]=curves21
  c_model_2[,,2]=curves22
  
  #3
  curves31=matrix(0,nrow = n,ncol = p)
  curves32=matrix(0,nrow = n,ncol = p)
  for (i in 1:n) {
    curves31[i,] =cos(40*pi*z)/10 +A[i,]
    curves32[i,] =sin(40*pi*z)/10 +B[i,]
  }
  c_model_3 = array(0,dim = c(n,p,2))
  c_model_3[,,1]=curves31
  c_model_3[,,2]=curves32
  
  #4 from Dai and Genton 2019
  #base4
  bcurves41=matrix(0,nrow = n,ncol = p)
  bcurves42=matrix(0,nrow = n,ncol = p)
  for (i in 1:n) {
    U_1=runif(1,2,3)
    U_2=runif(1,2,3)
    bcurves41[i,] =U_1*cos(80*pi*z)/10+A[i,]
    bcurves42[i,] =U_2*sin(40*pi*z)/10+B[i,]
  }
  bc_model_4 = array(0,dim = c(n,p,2))
  bc_model_4[,,1]=bcurves41
  bc_model_4[,,2]=bcurves42
  
  #contaminated model
  curves41=matrix(0,nrow = n,ncol = p)
  curves42=matrix(0,nrow = n,ncol = p)
  for (i in 1:n) {
    U_1=runif(1,3.2,3.5)
    U_2=runif(1,3.2,3.5)
    curves41[i,] =U_1*cos(4*pi*z)/10+A[i,]
    curves42[i,] =U_2*sin(4*pi*z)/10+B[i,]
  }
  c_model_4 = array(0,dim = c(n,p,2))
  c_model_4[,,1]=curves41
  c_model_4[,,2]=curves42
  
  #Create Models with outliers
  data_1=array(0,dim = c(n,p,2))
  data_2=array(0,dim = c(n,p,2))
  data_3=array(0,dim = c(n,p,2))
  data_4=array(0,dim = c(n,p,2))
  for (i in 1:n) {
    if(D[i]==1){
      data_1[i,,]=c_model_1[i,,]
      data_2[i,,]=c_model_2[i,,]
      data_3[i,,]=c_model_3[i,,]
      data_4[i,,]=c_model_4[i,,]
    }else{
      data_1[i,,]=data[i,,]
      data_2[i,,]=data[i,,]
      data_3[i,,]=data[i,,]
      data_4[i,,]=bc_model_4[i,,]
    }
    
  }
  
  
  mmodel=list(data_1,data_2, data_3,data_4)
  E=length(mmodel)
  
  #Outliers from MS plot
  out=list()
  for (m in 1:E) {
    msp= fdaoutlier::msplot(mmodel[[m]],plot=F)
    out[[m]]= msp$outliers
  }

  
  true_out=rep(0,E)
  false_out=rep(0,E)
  
  
  
  o = 1
  for (i in out) {
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
  FPR[,u]=(false_out/(n-total_out))*100
  
}


# True positive rate
mean.TPR=apply(TPR,1,mean);mean.TPR
sd.TPR=apply(TPR,1,sd);sd.TPR

# False positive rate
mean.FPR=apply(FPR,1,mean);mean.FPR
sd.FPR=apply(FPR,1,sd);sd.FPR
