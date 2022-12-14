#Code to repdroduce Figure 1. Magnitude outliers
library(fda)
set.seed(167)
M=100
N=100
k=6
#Generations of the models
cov.fun=function(d,k,c,mu){
  k*exp(-c*d^mu)
}
t=seq(0,1,len=M)
d=dist(t,upper=TRUE,diag=TRUE)
d.matrix=as.matrix(d)
C=rbinom(N,1,0.1)
#covariance function in time
t.cov=cov.fun(d.matrix,1,1,1)
# Cholesky Decomposition
L=chol(t.cov)
mu=0
e=matrix(rnorm(N*M),M,N)
#Model_0 Base Model
xt1=as.matrix(mu+t(L)%*%e)
#matplot(xt1[,C==1],type = 'l')
#Pure magnitude outliers
s=2*rbinom(N,1,0.5)-1
cs.m=matrix(C*s,M,N,byrow=TRUE)
xt2=matrix(NaN,nrow = M,ncol = N)
xt3=matrix(NaN,nrow = M,ncol = N)
xt4=matrix(NaN,nrow = M,ncol = N)
for (i in 1:N) {
  if(C[i]==1){
    xt2[,i]=as.matrix(xt1[,i]+k*cs.m[,i])
    
    #Partially Magnitude outliers
    xt3.=matrix(NaN,M,N)
    for (h in 1:N){
      TT=runif(1,0,1)
      xt3.[,h] = (t<TT) * xt1[,h] + (t>TT) * xt2[,h]
    }
    xt3[,i]=xt3.[,i]
    #Peaks magnitud outliers
    xt4.=matrix(NaN,M,N)
    for (h in 1:N) {
      l=0.5
      TTT=runif(1,0,1-l)
      xt4.[,h] =(t < TTT | t > (TTT + l)) * xt1[,h] + (t > TTT & t < (TTT+l)) * xt2[,h]
    }
    xt4[,i]=xt4.[,i]
  }else{
    xt2[,i]=xt1[,i]
    xt3[,i]=xt1[,i]
    xt4[,i]=xt1[,i]
  }
}

#Some plots of the magnitude-type outliers
par(mfrow=c(1,1))
matplot(xt1,type = "l",main="Functional data without outliers",ylab ="Simulated_data",xlab = "Time",col = 2,cex.main=0.8)
matplot(xt2,type = "l",main="Functional data with pure magnitud  outliers",ylab ="Simulated_data",xlab = "Time",col = 3,cex.main=0.7)
matplot(xt3,type = "l",main="Functional data with partially magnitud  outliers",ylab ="Simulated_data",xlab = "Time",col = 4,cex.main=0.7)
matplot(xt4,type = "l",main="Functional data with peaks magnitud  outliers",ylab ="Simulated_data",xlab = "Time",col=5,cex.main=0.7)
par(mfrow=c(1,1))

#Combination function
combinat=function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}
#########################################################################

# Check the functional boxplot of the magnitude outliers
model_1=list(xt1,xt2,xt3,xt4)
E=length(model_1)
out=list()
AllPWD=list()
mediancurve=list()
for (k in 1:E) {
  #Comparison with functional boxplot
  fda::fbplot(model_1[[k]],method = 'MBD',plot = F)
  out[[k]]= fda::fbplot(model_1[[k]],method = 'MBD',plot = F)$outpoint
  mediancurve[[k]]=which(fda::fbplot(model_1[[k]],method ='MBD',plot = F)$depth==max(fda::fbplot(model_1[[k]],method ='MBD',plot = F)$depth))
  #mediancurve2[[k]]=functional_boxplot(t(model_1[[k]]),depth_method = "mbd")$median_curve 
  #d[[k]]=functional_boxplot(model_1[[k]],depth_method = "mbd")$depth_values
  R1=as.matrix(t(apply(model_1[[k]], 1,rank))) 
  U=N-R1
  D=R1-1
  Point.wiseall=as.matrix(((U*D)+N-1)/combinat(N,2))
  AllPWD[[k]]=Point.wiseall
  #boxplot(AllPWD[[k]],col=k+1)
}
out
mediancurve
par(mfrow=c(1,1))
#Median Curves

#################################################################################
#Figure 1 of  the paper
#################################################################################


#Pure Magnitud
mixed_models=matrix(0,nrow = M,ncol = N)
for (i in 1:N) {
  if(i==48){
    mixed_models[,i]=xt2[,i]
  }else if(i==40){
    mixed_models[,i]=xt3[,i]
  }else if(i==19){
    mixed_models[,i]=xt4[,i]
  }else{
    mixed_models[,i]=xt1[,i]
  }
}
par(mfrow=c(1,1))
mod=mixed_models
med=mediancurve[[1]]
plot(NaN,xlim=c(0,M),ylim=c(min(mod),max(mod)),ylab='y',xlab='t')
for (i in 1:N) {
  lines(mod[,i],col="gray86")
}
lines(mod[,18],col=2,lwd=2)
lines(mod[,48],col=4,lwd=2)
lines(mod[,40],col=6,lwd=2)
lines(mod[,19],col=7)


mod1=AllPWD[[1]]
mod2=AllPWD[[2]]
mod3=AllPWD[[3]]
mod4=AllPWD[[4]]
boxplot(mod1[,18],mod2[,48],mod3[,40],mod4[,19],col=c(2,4,6,7),ylab='PWD')
#mtext("PWD of some magnitude outlier models", side = 3, line = -2, outer = TRUE)
par(mfrow=c(1,2))

