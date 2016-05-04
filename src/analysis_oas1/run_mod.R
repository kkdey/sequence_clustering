load("data/oas1/OAS1_mod.Robj")
source("src/cluster_seq_mix.R")

data = M

res = cluster.mix(data[, 2049:4096], smooth = TRUE, K = 4, tol = 1e-4, maxit = 4000)


