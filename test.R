source("cluster_seq.R")

pi.true=c(0.2,0.2,0.25,0.35)
lambda.true=matrix(0.1,nrow=4,ncol=1024)
lambda.true[1,100:200]=1
lambda.true[2,300:400]=1
lambda.true[3,500:600]=1
lambda.true[4,700:800]=1


y=matrix(0,nrow=20,ncol=1024)
clus.ind=0
for(i in 1:20){
  clus.ind[i]=sample(1:4,1,prob=pi.true)
  y[i,]=rpois(1024,lambda.true[clus.ind[i],])
}

clus.dis.res=cluster.dis(y,K=4,pseudocounts=1e-10,tol=1e-6,maxit=500)
clus.soft.res=cluster.soft(y,K=4,pseudocounts=1e-10,tol=1e-6,maxit=500)

kmeans.res=kmeans(y,4,nstart=5)

normalize(as.vector(table(kmeans.res$cluster)))
normalize(as.vector(table(clus.ind)))

library(smash)
##discrete
smooth.dis.res.pois=t(apply(clus.dis.res$lambda,1,ashsmooth.pois,cxx=FALSE))
smooth.dis.res.gaus=t(apply(clus.dis.res$lambda,1,ashsmooth.gaus))
y.gen=matrix(rpois(4*1024,clus.dis.res$lambda),nrow=4,ncol=1024)
smooth.dis.res.gen=t(apply(y.gen,1,ashsmooth.pois,cxx=FALSE))
plot(clus.dis.res$lambda[1,],pch=18,cex=0.5)
lines(smooth.dis.res.pois[1,],col=2)
lines(smooth.dis.res.gaus[1,],col=3)
lines(smooth.dis.res.gen[1,],col=4)


##soft
smooth.soft.res.pois=t(apply(clus.soft.res$lambda,1,ashsmooth.pois,cxx=FALSE))
smooth.soft.res.gaus=t(apply(clus.soft.res$lambda,1,ashsmooth.gaus))
y.gen=matrix(rpois(4*1024,clus.soft.res$lambda),nrow=4,ncol=1024)
smooth.soft.res.gen=t(apply(y.gen,1,ashsmooth.pois,cxx=FALSE))
plot(clus.soft.res$lambda[1,],pch=18,cex=0.5)
lines(smooth.soft.res.pois[1,],col=2)
lines(smooth.soft.res.gaus[1,],col=3)
lines(smooth.soft.res.gen[1,],col=4)


plot(clus.res$lambda[1,],pch=18,cex=0.5)

lines(smooth.res.pois[1,],col=4)



PoisMixClus(y,4,4,FALSE,conds=1:1024)
