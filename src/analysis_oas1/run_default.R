load("data/oas1/OAS1.Robj")
source("src/cluster_seq_mix.R")

data = Robj$M
data = as.matrix(data$M)


res = cluster.mix(data[, 2049:4096], smooth = TRUE, K = 4, tol = 1e-4, maxit = 4000)



