gene_region = read.table("supplemental/gtex/gene_region.txt", stringsAsFactors = FALSE)
gene_region = unlist(gene_region, use.names = FALSE)
gene_region = gsub("\\:", "_", gene_region)
gene_region = gsub("\\-", "_", gene_region)

#check to see if there are any regions with 0 total reads
for(i in 1:length(gene_region)){
  load(paste0("data/gtex/reads_100_", gene_region[i], ".Robj"))
  
  check = 0
  for(i in 1:8){
    check = check + sum(rowSums(reads[[i]][[2]]) == 0)
  }
  print(check)
}
