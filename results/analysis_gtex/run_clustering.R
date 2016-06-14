source("../../src/cluster_seq_mix.R")



args = commandArgs(TRUE)
K = as.numeric(args[1])
sample_size = as.numeric(args[2])
region = as.character(args[3])

if(sample_size !=0){
  file = paste0("reads_100_", region, ".Robj")
}else{
  file = paste0("reads_all_", region, ".Robj")
}
load(file.path("../../data/gtex", file))


data_matrix = NULL
tissue_name = NULL
for(i in 1:8){
  data_temp = reads[[i]][[2]]
  if(sample_size != 0){
    sample_ind = sample(1:100, sample_size)
    data_temp = data_temp[sample_ind, ]
  }
  data_matrix = rbind(data_matrix, data_temp)
  tissue_name_temp = reads[[i]][[1]]
  tissue_name_temp = gsub(" .*$", "", tissue_name_temp)
  tissue_name = c(tissue_name, tissue_name_temp)
}

res = cluster.mix(data_matrix, smooth = TRUE, K = K, tol = 1e-4, maxit = 10000)

if(sample_size != 0){
  save_name = paste("cluster_res_100", K, "clusters", region, sep = "_")
}else{
  save_name = paste("cluster_res_all", K, "clusters", region, sep = "_")
}
save_name = paste0(save_name, ".RData")
save.image(save_name)