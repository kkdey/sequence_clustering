normalize=function(x) x/sum(x)

smooth.pois=function(yr){
  lambda.smooth=t(apply(yr,1,ashsmooth.pois))
}



####discrete mixture
EMupd.dis=function(y,pi,lambda,n,K,B,pseudocounts){
  #gamma is n*K, pi is 1*K, lambda is K*B, y is n*B
  lik.ini=matrix(1,nrow=n,ncol=B)%*%t(-lambda)+y%*%log(t(lambda))-log(factorial(y))%*%matrix(1,nrow=B,ncol=K)
  lik.ini=lik.ini-max(lik.ini)  #subtract row max?
  gamma=(rep(1,n)%o%pi)*exp(lik.ini)
  gamma=t(apply(gamma,1,normalize))
  pi=colSums(gamma)/n
  lambda=(t(gamma)%*%y)/(t(gamma)%*%matrix(1,nrow=n,ncol=B))
  lambda[is.na(lambda)]=0
  lambda[lambda==0]=pseudocounts  #any way to avoid this by having 0*log(0)=1?
  return(list(pi=pi,lambda=lambda,gamma=gamma))
}

negloglik.dis=function(y,pi,lambda,n,K,B){
  loglik.ini=matrix(1,nrow=n,ncol=B)%*%t(-lambda)+y%*%log(t(lambda))-log(factorial(y))%*%matrix(1,nrow=B,ncol=K)
  loglik.ini.max=max(loglik.ini)
  loglik.ini=loglik.ini-loglik.ini.max
  loglik.ind=exp(loglik.ini)
  loglik.tot=-sum(loglik.ini.max+log(loglik.ind%*%pi))
  return(loglik.tot)
}



EMproc.dis=function(y,pi,lambda,n,K,B,pseudocounts,tol,maxit){
  loglik.old=Inf
  loglik=negloglik.dis(y,pi,lambda,n,K,B)
  cyc=0
  while(abs(loglik-loglik.old)>tol&cyc<maxit){
    loglik.old=loglik
    res=EMupd.dis(y,pi,lambda,n,K,B,pseudocounts)
    pi=res$pi
    lambda=res$lambda
    gamma=res$gamma
    loglik=negloglik.dis(y,pi,lambda,n,K,B)
    cyc=cyc+1
#print(cyc)
#print(pi)
#print(abs(loglik-loglik.old))
  }
  return(list(pi=pi,lambda=lambda,gamma=gamma,loglik=loglik))
}

cluster.dis=function(y,pi0=NULL,lambda0=NULL,K,pseudocounts,tol,maxit){
  n=dim(y)[1]
  B=dim(y)[2]
  
  if(is.null(pi0)|is.null(lambda0)){
    kmeans.init=kmeans(y,K,nstart=5)
  }

  if(is.null(pi0)) pi0=normalize(as.vector(table(kmeans.init$cluster)))

  if(is.null(lambda0)){
    lambda0=kmeans.init$centers
    row.names(lambda0)=NULL
  }
  lambda0[lambda0==0]=pseudocounts 

  out=EMproc.dis(y,pi0,lambda0,n,K,B,pseudocounts,tol,maxit)
  return(list(pi=out$pi,lambda=out$lambda,gamma=out$gamma,loglik=out$loglik))
}



####soft clustering
EMupd.soft=function(y,pi,lambda,n,K,B){
  #gamma is nB*K, pi is n*K, lambda is K*B, y is n*B
  #_.ext are all nB*K
  pi.ext=pi[rep(1:n,each=B),]
  lambda.ext=t(lambda)[rep(1:B,n),]
  y.ext=as.vector(t(y))%o%rep(1,K)
  yll=y.ext*log(lambda.ext)
  yll[is.na(yll)]=0
  gamma.ini=log(pi.ext)-lambda.ext+yll-log(factorial(y.ext))
  gamma=exp(gamma.ini-max(gamma.ini)) #subtract row max?
  gamma=t(apply(gamma,1,normalize))
  pi=(diag(1,n)[,rep(1:n,each=B)]%*%gamma)/B
  lambda=t((diag(1,B)[,rep(1:B,n)]%*%(gamma*y.ext))/(diag(1,B)[,rep(1:B,n)]%*%gamma))
  lambda[is.na(lambda)]=0
  return(list(pi=pi,lambda=lambda,gamma=gamma))
}

negloglik.soft=function(y,pi,lambda,n,K,B){
  pi.ext=pi[rep(1:n,each=B),]
  lambda.ext=t(lambda)[rep(1:B,n),]
  y.ext=as.vector(t(y))%o%rep(1,K)
  yll=y.ext*log(lambda.ext)
  yll[is.na(yll)]=0
  loglik.ini=log(pi.ext)-lambda.ext+yll-log(factorial(y.ext))
  loglik.ini.max=max(loglik.ini)
  loglik.ini=loglik.ini-loglik.ini.max
  loglik.ind=exp(loglik.ini)
  loglik.tot=-sum(loglik.ini.max+log(rowSums(loglik.ind)))
  return(loglik.tot)
}



EMproc.soft=function(y,pi,lambda,n,K,B,tol,maxit){
  loglik.old=Inf
  loglik=negloglik.soft(y,pi,lambda,n,K,B)
  cyc=0
  while(abs(loglik-loglik.old)>tol&cyc<maxit){
    loglik.old=loglik
    res=EMupd.soft(y,pi,lambda,n,K,B)
    pi=res$pi
    lambda=res$lambda
    gamma=res$gamma
    loglik=negloglik.soft(y,pi,lambda,n,K,B)
    cyc=cyc+1
#print(cyc)
#print(pi)
print(loglik)
  }
  return(list(pi=pi,lambda=lambda,gamma=gamma,loglik=loglik))
}

cluster.soft=function(y,pi0=NULL,lambda0=NULL,K,pseudocounts,tol,maxit){
  n=dim(y)[1]
  B=dim(y)[2]
  
  if(is.null(pi0)|is.null(lambda0)){
    kmeans.init=kmeans(y,K,nstart=5)
  }

  if(is.null(pi0)) pi0=rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))

  if(is.null(lambda0)){
    lambda0=kmeans.init$centers
    row.names(lambda0)=NULL
  }
  lambda0[lambda0==0]=pseudocounts 

  out=EMproc.soft(y,pi0,lambda0,n,K,B,tol,maxit)
  return(list(pi=out$pi,lambda=out$lambda,gamma=out$gamma,loglik=out$loglik))
}


# EMupd.soft.check=function(y,pi,lambda,n,K,B){
#   gamma=array(0,dim=c(n,K,B))
#   for(i in 1:n){
#     for(k in 1:K){
#       for(b in 1:B){
#         gamma[i,k,b]=log(pi[i,k])-lambda[k,b]+y[i,b]*log(lambda[k,b])-log(factorial(y[i,b]))
#       }
#     }
#   }
#   
#   gamma=exp(gamma-max(gamma)) #subtract row max?
#   gamma=aperm((apply(gamma,c(1,3),normalize)),c(2,1,3))
#   pi=apply(gamma,c(1,2),sum)/B
#   gy=gamma
#   for(i in 1:n){
#     for(k in 1:K){
#       for(b in 1:B){
#         gy[i,k,b]=gamma[i,k,b]*y[i,b]
#       }
#     }
#   }
#   lambda=apply(gy,c(2,3),sum)/apply(gamma,c(2,3),sum)
#   lambda[is.na(lambda)]=0
#   return(list(pi=pi,lambda=lambda,gamma=gamma))
# }

# negloglik.soft.check=function(y,pi,lambda,n,K,B){
#   loglik.ini=array(0,dim=c(n,K,B))
#   
#   for(i in 1:n){
#     for(k in 1:K){
#       for(b in 1:B){
#         loglik.ini[i,k,b]=log(pi[i,k])-lambda[k,b]+y[i,b]*log(lambda[k,b])-log(factorial(y[i,b]))
#       }
#     }
#   }
#   loglik.ini.max=max(loglik.ini)
#   loglik.ini=loglik.ini-loglik.ini.max
#   loglik.ind=exp(loglik.ini)
#   loglik.tot=-sum(loglik.ini.max+log(apply(loglik.ind,c(1,3),sum)))
#   return(loglik.tot)
# }

