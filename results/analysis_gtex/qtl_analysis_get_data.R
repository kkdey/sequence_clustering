path = "../../data/gtex"
tissue_snp_file_names = list.files(path, ".txt.gz")
tissue_name_short = gsub("_.*$", "", tissue_snp_file_names)

gene_region = read.table("../../supplemental/gtex/gene_region.txt", stringsAsFactors = FALSE)
gene_region = unlist(gene_region, use.names = FALSE)
gene_region = gsub("\\:", "_", gene_region)
gene_region = gsub("\\-", "_", gene_region)
gene_region = strsplit(gene_region, "_")

j = 4
radius = 100000


for(k in 1:length(tissue_snp_file_names)){
  
  gene_chr = gene_region[[j]][1]
  gene_start = gene_region[[j]][2]
  gene_end = gene_region[[j]][3]
  range_start = as.numeric(gene_start) - radius
  range_end = as.numeric(gene_end) + radius
  
  con = file(file.path(path, tissue_snp_file_names[k]), open = "r")
  
  oneLine <- readLines(con, n = 1, warn = FALSE)
  line = strsplit(oneLine, "\t")
  line = line[[1]]
  subject_id = line[-1]
  
  geno = NULL
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0){
     line = strsplit(oneLine, "\t")
     line = line[[1]]
     snp_id = line[1]
     snp_id = strsplit(snp_id, "_")
     snp_chr = snp_id[[1]][1]
     snp_loc = as.numeric(snp_id[[1]][2])
     if(snp_chr != gene_chr | snp_loc < range_start | snp_loc > range_end) next
     geno = rbind(geno, line[-1])
     rownames(geno)[dim(geno)[1]] = line[1]
  }
  
  close(con)
  
  colnames(geno) = subject_id
  
  geno = as.data.frame(geno)
  geno = cbind(id = rownames(geno), geno)
  
#   write(t(geno), file = file.path(path, paste0(paste("genotype", paste(gene_region[[j]], collapse = "_"), tissue_name_short[k], sep = "_"), ".txt")), ncolumns = length(subject_id))
  save(geno, file = file.path(path, paste0(paste("genotype", paste(gene_region[[j]], collapse = "_"), tissue_name_short[k], sep = "_"), ".Robj")))
}