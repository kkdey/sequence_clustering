gene_region = read.table("../../supplemental/gtex/gene_region.txt", stringsAsFactors = FALSE)
gene_region = unlist(gene_region, use.names = FALSE)
gene_region = gsub("\\:", "_", gene_region)
gene_region = gsub("\\-", "_", gene_region)


K = 2

for(i in 1:length(gene_region)){
  load_name = paste("cluster_res_100", K, "clusters", gene_region[i], sep = "_")
  load_name = paste0(load_name, ".RData")
  load(load_name)
  
  n = dim(res$pi)[1]
  sep_lines = seq(0, n, length.out = 9)
  sep_lines = sep_lines[c(-1, -length(sep_lines))]
  
  par(mfrow = c(2, 1))
  plot(res$phi[1, ], type = 'l', col = 2)
  for(k in 2:K){
    lines(res$phi[k, ], col = k + 1)
  }
  barplot(t(res$pi), col = 2:3, axisnames = F, space = 0, border = NA, main = "structure plot for smoothed model", las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
  abline(v = sep_lines, lwd = 2)
}
