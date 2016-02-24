library(smash)
load("data/oas1/OAS1.Robj")
source("src/cluster_seq_mix.R")


data = Robj$M
data = as.matrix(data$M)


data_smooth = apply(data, 1, ashsmooth.pois) 
data_smooth = t(data_smooth)

res = cluster.mix(data_smooth[, 2049:4096], smooth = FALSE, K = 3, tol = 1e-4, maxit = 4000)


