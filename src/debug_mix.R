kmeans.init=kmeans(y,K,nstart=5)
pi0=rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
lambda0=kmeans.init$centers
lambda0[lambda0==0]=1e-6
phi0=t(apply(lambda0,1,normalize))
negloglik.mix(y,pi0,phi0,n,K,B)
tt=EMupd.mix(y,pi0,phi0,n,K,B)

gamma=pi0[rep(1:n,each=B),]*t(phi0)[rep(1:B,n),]
  gamma=t(apply(gamma,1,normalize))
  gammab=(as.vector(t(y))%o%rep(1,K))*gamma
  pi.num=(diag(1,n)[,rep(1:n,each=B)])%*%gammab
  pi=pi.num/(rowSums(y)%o%rep(1,K))
  ybt=(diag(1,B)[,rep(1:B,n)])%*%gammab
  phi=t(ybt/(rep(1,B)%o%colSums(gammab)))

gamma=res$gamma
gammab=(as.vector(t(y))%o%rep(1,K))*gamma
  ybt=(diag(1,B)[,rep(1:B,n)])%*%gammab
  ybw=(diag(1,B)[,rep(1:B,n)])%*%gamma
  ykb=ybt/ybw
  ykb[is.na(ykb)]=0
  ykt=colSums(ykb)

colSums(ybt)/colSums(res$pi)

colSums(ybt)[3]/sum(pi.true[,2])