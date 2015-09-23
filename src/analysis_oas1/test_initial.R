load("data/oas1/OAS1.Robj")
source("src/cluster_seq_mix.R")

data = Robj$M
data = as.matrix(data$M)


res=cluster.mix(data[, 2049:4096],smooth=TRUE,K=3,tol=1e-4,maxit=4000)


plot(res$phi[2,],type='l')
lines(res$phi[1,],col=2)
lines(res$phi[3,],col=4)


plot(res$lambda[2,],type='l')
lines(res$lambda[1,],col=2)
lines(res$lambda[3,],col=4)
