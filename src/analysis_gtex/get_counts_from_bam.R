
library(multiseq)

setwd("~/projects/sequence_clustering")
source("src/analysis_gtex/get_counts_single.R")


args = commandArgs(TRUE)
sample_size = as.numeric(args[1])
region = as.character(args[2])


region_split = split_region(region)
region_split$chr = gsub("[^0-9]", "", region_split$chr)

load("supplemental/gtex/runinfo_subset.Robj")


#sample_fixed = function(x) sample(x, size = sample_size)

samples_subset = aggregate(Run_s ~ body_site_s, data = runinfo_subset, FUN = sample)



reads = list()

for(i in 1:8){
  reads[[i]] = list()
  reads[[i]][[1]] = samples_subset[i, ][[1]]
#   reads[[i]][[2]] = matrix(0, nr = sample_size, nc = region_split$end - region_split$start + 1)
  temp = NULL
  j = 1
  k = 1
  if(sample_size != 0){
    while(j <= sample_size){
      filename = paste(samples_subset[i, ][[2]][[1]][k], region_split$chr, region_split$start, region_split$end, sep = "_")
      filename = paste0(filename, ".bam")
      bamfile = file.path("/mnt/gluster/data/external_private_supp/ncbi/dbGaP-3253", filename)
      if(file.exists(bamfile)){
        counts = get.counts.single(bamfile, region)
  #       reads[[i]][[2]][j, ] = counts
        temp = rbind(temp, counts)       
        if(sum(counts) == 0){
          k = k + 1
          next
        }
        j = j + 1
        k = k + 1
      }else{
        k = k + 1
        next
      }
    } 
  }else{
    while(k <= length(samples_subset[1, ][[2]][[1]])){
      filename = paste(samples_subset[i, ][[2]][[1]][k], region_split$chr, region_split$start, region_split$end, sep = "_")
      filename = paste0(filename, ".bam")
      bamfile = file.path("/mnt/gluster/data/external_private_supp/ncbi/dbGaP-3253", filename)
      if(file.exists(bamfile)){
        counts = get.counts.single(bamfile, region)
        #       reads[[i]][[2]][j, ] = counts
        if(sum(counts) == 0){
          k = k + 1
          next
        }
        temp = rbind(temp, counts)       
        k = k + 1
      }else{
        k = k + 1
        next
      }
    }
  }
  reads[[i]][[2]] = temp
}

reads[[9]] = region

if(sample_size != 0){
  save_name = paste("reads", sample_size, region_split$chr, region_split$start, region_split$end, sep = "_")
}else{
  save_name = paste("reads_all", region_split$chr, region_split$start, region_split$end, sep = "_")
}
save_name = paste0(save_name, ".Robj")
save(reads, file = file.path("data/gtex", save_name))