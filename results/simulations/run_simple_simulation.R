
args = commandArgs(TRUE)
i = as.integer(args[1])

mse = function(x, y) mean((x - y)^2)

source("../../src/cluster_seq_mix.R")

get.max=function(x) which(x==max(x))
normalize=function(x){
  #if(sum(abs(x))!=0){
  return(x/sum(x))
  #}else{
  #  return(rep(0,length(x)))
  #}
}


#generate data

K=4
n=20
B=1024
lambda.true=matrix(0.1,nrow=4,ncol=1024)
lambda.true[1,100:200]=1
lambda.true[2,300:400]=10
lambda.true[3,500:600]=3  #ML doesn't work well if this is 1
lambda.true[4,700:800]=5
phi.true=t(apply(lambda.true,1,normalize))


  set.seed(10*i)
  
  pi.true=matrix(0,n,4)
  for(j in 1:n){
    clus.ini=sample(1:K,1)
    pi.true[j,clus.ini]=0.7
    pi.true[j,-clus.ini]=0.1
  }
  
  y=matrix(0,nrow=n,ncol=B)
  for(k in 1:n){
    for(b in 1:B){
      lambda=sum(pi.true[k,]*lambda.true[,b])
      y[k,b]=rpois(1,lambda)
    }
  }
  
  
  #run different methods
  
  res.nosmooth=cluster.mix(y,smooth=FALSE,K=4,tol=1e-4,maxit=500)
  res.smooth=cluster.mix(y,smooth=TRUE,K=4,tol=1e-4,maxit=800)
  
  mse.pi.smooth=mse(pi.true,res.smooth$pi)
  mse.phi.smooth=mse(phi.true,res.smooth$phi)
  mse.pi.nosmooth=mse(pi.true,res.nosmooth$pi)
  mse.phi.nosmooth=mse(phi.true,res.nosmooth$phi)
  

save.image(paste0("simple_simulation_mse_", i, ".RData"))