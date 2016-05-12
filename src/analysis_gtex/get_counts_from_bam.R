
library(multiseq)
source("src/analysis_gtex/get_counts_single.R")


args = commandArgs(TRUE)
sample_size = as.numeric(args[1])
region = as.character(args[2])


region_split = split_region(region)
region_split$chr = gsub("[^0-9]", "", region_split$chr)

load("supplemental/gtex/runinfo_subset.Robj")


sample_fixed = function(x) sample(x, size = sample_size)

samples_subset = aggregate(Run_s ~ body_site_s, data = runinfo_subset, FUN = sample_fixed)



reads = list()

for(i in 1:8){
  reads[[i]][[1]] = samples_subset[i, ][[1]]
  reads[[i]][[2]] = matrix(0, nr = sample_size, nc = region_split$end - region_split$start + 1)
  for(j in 1:sample_size){
    filename = paste(samples_subset[i, ][[2]][1, j], region_split$chr, region_split$start, region_split$end, sep = "_")
    filename = paste0(filename, ".bam")
    counts = get_counts_single(filename, region)
    reads[[i]][[2]][j, ] = counts    
  } 
}

reads[[9]] = region

save_name = paste("reads", sample_size, region_split$chr, region_split$start, region_split$end, sep = "_")
save_name = paste0(save_name, ".Robj")
save(reads, file = file.path("data/gtex", save_name))