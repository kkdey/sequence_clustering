mse.pi.smooth.all=0
mse.phi.smooth.all=0
mse.pi.nosmooth.all=0
mse.phi.nosmooth.all=0

for(i in 1:100){
  load(paste0("simple_simulation_mse_", i, ".RData"))
  ordering.smooth=NULL
  ordering.nosmooth=NULL
  for(k in 1:4){
    mse.temp.smooth=apply((res.smooth$phi-rep(1,4)%o%phi.true[k,])^2,1,mean)
    mse.temp.nosmooth=apply((res.nosmooth$phi-rep(1,4)%o%phi.true[k,])^2,1,mean)
    ordering.smooth[k]=which(mse.temp.smooth==min(mse.temp.smooth))
    ordering.nosmooth[k]=which(mse.temp.nosmooth==min(mse.temp.nosmooth))    
  }
  
  mse.pi.smooth.all[i]=mse(pi.true,res.smooth$pi[,ordering.smooth])
  mse.phi.smooth.all[i]=mse(phi.true,res.smooth$phi[ordering.smooth,])
  mse.pi.nosmooth.all[i]=mse(pi.true,res.nosmooth$pi[,ordering.nosmooth])
  mse.phi.nosmooth.all[i]=mse(phi.true,res.nosmooth$phi[ordering.nosmooth,])
}

summary(mse.pi.smooth.all)
summary(mse.pi.nosmooth.all)
summary(mse.phi.smooth.all)
summary(mse.phi.nosmooth.all)

sd(mse.pi.smooth.all)
sd(mse.pi.nosmooth.all)
sd(mse.phi.smooth.all)
sd(mse.phi.nosmooth.all)

pdf('clustering_simulation_example.pdf',height=10,width=8)
par(mfrow=c(4,1),mar=c(5,4,1,2),oma=c(2,2,0.4,0.4))
plot(res.nosmooth$phi[ordering.nosmooth[1],],xlab='location',ylab='normalized intensity',type='l',col=3)
lines(phi.true[1,],col=1,lwd=2)
lines(res.smooth$phi[ordering.smooth[1],],col=2,lwd=2)
legend("topright",legend=c("True profile", "Cluster-seq", "Basic GoM model"), lty=1, col=1:3)
plot(res.nosmooth$phi[ordering.nosmooth[2],],xlab='location',ylab='normalized intensity',type='l',col=3)
lines(phi.true[2,],col=1,lwd=2)
lines(res.smooth$phi[ordering.smooth[2],],col=2,lwd=2)
legend("topright",legend=c("True profile", "Cluster-seq", "Basic GoM model"), lty=1, col=1:3)
plot(res.nosmooth$phi[ordering.nosmooth[3],],xlab='location',ylab='normalized intensity',type='l',col=3)
lines(phi.true[3,],col=1,lwd=2)
lines(res.smooth$phi[ordering.smooth[3],],col=2,lwd=2)
legend("topright",legend=c("True profile", "Cluster-seq", "Basic GoM model"), lty=1, col=1:3)
plot(res.nosmooth$phi[ordering.nosmooth[4],],xlab='location',ylab='normalized intensity',type='l',col=3)
lines(phi.true[4,],col=1,lwd=2)
lines(res.smooth$phi[ordering.smooth[4],],col=2,lwd=2)
dev.off()
