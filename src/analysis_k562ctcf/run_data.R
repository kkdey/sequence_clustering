dir.name = "~/projects/sequence_clustering"
source(file.path(dir.name,"src","cluster_seq_mix.R"))
library(smash)

load(file.path(dir.name,"results/analysis_k562ctcf","data.sig.Robj"))

ind.high = rowMeans(data.sig)>mean(rowMeans(data.sig))
sample.data = data.sig[ind.high, ][sample(1:sum(ind.high), 2000), ]
print(dim(sample.data))
sample.res.smooth = cluster.mix(sample.data, smooth=TRUE, K=10, tol=1e-4, maxit=10000)
