dir.name = "~/projects/sequence_clustering"
source(file.path(dir.name,"src","cluster_seq_mix.R"))
library(smash)

load(file.path(dir.name,"results/analysis_k562ctcf","data.sig.Robj"))

sample.data = data.sig[sample(1:dim(data.sig)[1], 100), ]
sample.res.smooth = cluster.mix(sample.data, smooth=TRUE, K=10, tol=1e-4, maxit=10000)
