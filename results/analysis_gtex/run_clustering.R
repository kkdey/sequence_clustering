source("../../cluster_seq_mix.R")



args = commandArgs(TRUE)
K = args[1]
region = args[2]

file = paste0("reads_100_", region, ".Robj")
load(file.path("../../data/gtex", file))


data_matrix = NULL
tissue_name = NULL
for(i in 1:8){
  data_temp = reads[[i]][[2]]
  sample_ind = sample(1:100, 5)
  data_temp = data_temp[sample_ind, ]
  data_matrix = rbind(data_matrix, data_temp)
  tissue_name_temp = reads[[i]][[1]]
  tissue_name_temp = gsub(" .*$", "", tissue_name_temp)
  tissue_name = c(tissue_name, tissue_name_temp)
}

res = cluster.mix(data_matrix, smooth = TRUE, K = K, tol = 1e-4, maxit = 5000)

save_name = paste("cluster_res_100", K, "clusters", region)
save_name = paste0(save_name, ".RData")
save.image(save_name)