dir.name = "~/projects/sequence_clustering"
source(file.path("src","cluster_seq_mix.R"))
library(smash)

load(file.path(dir.name,"results/analysis_k562ctcf","data.sig.Robj"))


ind.high = rowMeans(data.sig)>mean(rowMeans(data.sig))
sample.data = data.sig[ind.high, ][sample(1:sum(ind.high), 500), ]
sample.data = data.sig[ind.high, ]
sample.res.smooth = cluster.mix(sample.data, smooth=TRUE, K=10, tol=1e-4, maxit=10000)


for(i in 1:10){
  plot(sample.res.smooth$lambda[i,],type='l',main=paste0("cluster ",i))
}


sample.res.smooth$pi
