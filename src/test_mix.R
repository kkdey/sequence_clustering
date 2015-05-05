library(maptpx)
library(smash)
get.max=function(x) which(x==max(x))
normalize=function(x){
  #if(sum(abs(x))!=0){
    return(x/sum(x))
  #}else{
  #  return(rep(0,length(x)))
  #}
}


K=4
n=20
B=1024
lambda.true=matrix(0.1,nrow=4,ncol=1024)
lambda.true[1,100:200]=1
lambda.true[2,300:400]=10
lambda.true[3,500:600]=3  #ML doesn't work well if this is 1
lambda.true[4,700:800]=5
phi.true=t( apply(lambda.true,1,normalize))

pi.true=matrix(0,n,4)
for(i in 1:n){
  set.seed(10*i)
  clus.ini=sample(1:K,1)
  pi.true[i,clus.ini]=0.7
  pi.true[i,-clus.ini]=0.1
}


yt=0
clus.ind=0
y=matrix(0,nrow=n,ncol=B)
for(i in 1:n){
  yt[i]=rpois(1,sum(pi.true[i,]*rowSums(lambda.true)))
  p.mult=colSums((pi.true[i,]%o%rep(1,B))*phi.true)
  y[i,]=rmultinom(1,yt[i],p.mult)
}


y=matrix(0,nrow=n,ncol=B)
for(i in 1:n){
  for(b in 1:B){
    lambda=sum(pi.true[i,]*lambda.true[,b])
    y[i,b]=rpois(1,lambda)
  }
}

res=cluster.mix(y,smooth=FALSE,K=4,tol=1e-4,maxit=4000)
res.smooth=cluster.mix(y,K=4,tol=1e-4,maxit=4000)
res.tpx=topics(y,4)
apply(res.tpx$omega,1,get.max)
apply(res$pi,1,get.max)

par(mfrow=c(4,1))
plot(res$phi[1,],type='l')
lines(res.tpx$theta[,2],col=2)
lines(res.smooth$phi[1,],col=4)
plot(res$phi[2,],type='l')
lines(res.tpx$theta[,1],col=2)
lines(res.smooth$phi[2,],col=4)
plot(res$phi[3,],type='l')
lines(res.tpx$theta[,3],col=2)
lines(res.smooth$phi[3,],col=4)
plot(res$phi[4,],type='l')
lines(res.tpx$theta[,4],col=2)
lines(res.smooth$phi[4,],col=4)





lambda.smooth=t(apply(res$lambda,1,ashsmooth.pois,cxx=FALSE))
par(mfrow=c(4,1))
plot(lambda.smooth[1,],ylim=c(0,10),type='l',col=4)
lines(res.smooth$lambda[1,],col=2)
lines(lambda.true[2,])
plot(lambda.smooth[2,],ylim=c(0,5),type='l',col=4)
lines(res.smooth$lambda[2,],col=2)
lines(lambda.true[4,])
plot(lambda.smooth[3,],ylim=c(0,4),type='l',col=4)
lines(res.smooth$lambda[4,],col=2)
lines(lambda.true[3,])
plot(lambda.smooth[4,],ylim=c(0,2),type='l',col=4)
lines(res.smooth$lambda[3,],col=2)
lines(lambda.true[1,])
