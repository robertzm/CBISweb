#########################################################
#
#This R function provides the nonparametric estimator of Lin's CCC 
#  to run the function:  nonpar.CCC.est(data)
# input: data=data.frame(cbind(X1,X2,eta1,eta2))
#  X1,X2  observed survival times based on two raters/methods 
#  eta1,eta2, censoring indicator which equals 1 is X1,X2 are observed event times and 0 if X1 and X2 
#  represents censoring time.
#  
#########################################################


#functions for calculating the true values

SC.func=function(T,censor){
deathtime=unique(sort(T[censor[]==1]))
ndeath=as.numeric(table(sort(T[censor[]==1])))
nrisk=rep(0,length(deathtime))
for(i in 1:length(deathtime)){nrisk[i]=sum((deathtime[i]<=T))}
prodobj=1-ndeath/nrisk
survp=rep(0,length(deathtime))
for(i in 1:length(deathtime)){survp[i]=prod(prodobj[1:i])}
return(data.frame(cbind(deathtime,ndeath,nrisk,survp)))}

biv.survival.ly.nosd=function(da,t0,SC){
#t0 is a m*2 matrix of survival times
#SC=stepfun(SC0$deathtime,c(1,SC0$survp),right=TRUE)
m=dim(t0)[1]
n=dim(da)[1]
inv.prob=SC(apply(t0,1,max))
X1.mat=matrix(rep(da$X1,m),n)
X2.mat=matrix(rep(da$X2,m),n)
t0.1.mat=matrix(rep(t0[,1],rep(n,m)),n)
t0.2.mat=matrix(rep(t0[,2],rep(n,m)),n)
compare.mat=(X1.mat>t0.1.mat)&(X2.mat>t0.2.mat)
mean.mat=apply(compare.mat,2,mean)
surv.prob=mean.mat/inv.prob
return(list(surv=surv.prob))}

fisherz=function(rho){
fisherz=0.5*log((1+rho)/(1-rho))
return(fisherz)}

dfisherz=function(rho){
dfisherz=1/(1-rho^2)
return(dfisherz)}

nonpar.CCC.est=function(da){
C.time=pmax(da$X1,da$X2)
C.censor=1-da$eta1*da$eta2
SC0=SC.func(C.time,C.censor)
if(sum(C.censor)>0){SC=stepfun(SC0$deathtime,c(1,SC0$survp),right=TRUE)}
if(sum(C.censor)==0){SC=stepfun(c(0,100),c(1,1,1),right=TRUE)}

Cy=function(t){Cy=ecdf(rep(0,n))(t)-ecdf(C.time)(t-1e-10)}
max.tau=max(da$X1,da$X2)+0.001

#first the marginals
t1.times=c(0,sort(unique(da$X1)),max.tau)
t2.times=c(0,sort(unique(da$X2)),max.tau)
lt1=length(t1.times)
lt2=length(t2.times)
t1.interval=t1.times[-1]-t1.times[-lt1]
t2.interval=t2.times[-1]-t2.times[-lt2]
surv.1=biv.survival.ly.nosd(da,cbind(t1.times,0),SC)
surv.2=biv.survival.ly.nosd(da,cbind(0,t2.times),SC)

surv.1.adj=surv.1$surv[-lt1]
surv.2.adj=surv.2$surv[-lt2]
t1.adj=t1.times[-lt1]
t2.adj=t2.times[-lt2]

ET1=sum(surv.1.adj*t1.interval)
ET2=sum(surv.2.adj*t2.interval)
EtT1=sum(surv.1.adj*t1.adj*t1.interval)
EtT2=sum(surv.2.adj*t2.adj*t2.interval)
var.T1=2*EtT1-ET1^2
var.T2=2*EtT2-ET2^2

#construct a grid based on observed survival times
t.grid=t1.grid=t.area=c()
row.space=t1.times[2:length(t1.times)]-t1.times[1:(length(t1.times)-1)]
col.space=t2.times[2:length(t2.times)]-t2.times[1:(length(t2.times)-1)]

for(i in 1:(length(t1.times)-1)){
  t.grid=rbind(t.grid,cbind(t1.times[i],t2.times[-length(t2.times)]))
  t.area=c(t.area,row.space[i]*col.space)
  }
  
surv.grid=biv.survival.ly.nosd(da,t.grid,SC)
cov12=sum(surv.grid$surv*t.area)-ET1*ET2
rho.est=(2*cov12)/(var.T1+var.T2+(ET1-ET2)^2)
return(rho.est)}




