negloglik.dis(y,rep(1/4,4),rep(1,4)%o%colMeans(y),50,4,1024)

lambda=rep(1,4)%o%(colMeans(y))
lambda[lambda==0]=1e-6
lambda=lambda+runif(K*B,0,0.001)
pi=normalize(runif(4))
n=50
K=4
B=1024
gamma=(rep(1,n)%o%pi)*exp(matrix(1,nrow=n,ncol=B)%*%t(-lambda)+y%*%log(t(lambda))-log(factorial(y))%*%matrix(1,nrow=B,ncol=K))
gamma=t(apply(gamma,1,normalize))
colSums(gamma)/n

y[1,]*log(t(lambda))[,1]
y[1,]*log(t(lambda))[,2]

negloglik.dis(y,pi.true,lambda.true,50,4,1024)
EMupd.dis(y,pi.true,lambda.true,n,K,B,1e-6)

gamma=(rep(1,n)%o%pi.true)*exp(matrix(1,nrow=n,ncol=B)%*%t(-lambda.true)+y%*%log(t(lambda.true))-log(factorial(y))%*%matrix(1,nrow=B,ncol=K))
gamma=t(apply(gamma,1,normalize))
colSums(gamma)/n


tt=kmeans(y,4,nstart=5)
colSums(y[tt$cluster==1,])/sum(tt$cluster==1)-tt$center[1,]
tt$center[1,]

negloglik.dis(y,normalize(as.vector(table(kmeans.res$cluster))),kmeans.res$centers,20,4,1024)