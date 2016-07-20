library(forecast)

gene_region = read.table("../../supplemental/gtex/gene_region.txt", stringsAsFactors = FALSE)
gene_region = unlist(gene_region, use.names = FALSE)
gene_region = gsub("\\:", "_", gene_region)
gene_region = gsub("\\-", "_", gene_region)


K = c(3, 4, 3)

gi = c(17, 2, 4)

g = 3


load_name = paste("reads_all", gene_region[gi[g]], sep = "_")
load_name = paste0(load_name, ".Robj")
load(file.path("../../data/gtex", load_name))
load("../../supplemental/gtex/runinfo_subset.Robj")


normalized_reads = list()
order_reads = list()
for(j in 1:(length(reads) - 1)){
  library_size = runinfo_subset[runinfo_subset$Run_s %in% run[[j]], ]$MBases_l*10^6/76
  normalized_reads[[j]] = rowSums(reads[[j]][[2]])/library_size
  order_reads[[j]] = order(normalized_reads[[j]])
}


gene_region = read.table("../../supplemental/gtex/gene_region.txt", stringsAsFactors = FALSE)
gene_region = unlist(gene_region, use.names = FALSE)
gene_region = gsub("\\:", "_", gene_region)
gene_region = gsub("\\-", "_", gene_region)

  load_name = paste("cluster_res_all", K[g], "clusters", gene_region[gi[g]], sep = "_")
  load_name = paste0(load_name, ".RData")
  if(!file.exists(load_name)) next
  load(load_name)
  
  sep_lines = 0
  for(j in 1:(length(reads) - 1)){
    temp = dim(reads[[j]][[2]])[1]
    sep_lines[j + 1] = sep_lines[j] + temp
  }
  
  sep_lines_mid = (sep_lines[1:(length(sep_lines) - 1)] + sep_lines[2:length(sep_lines)])/2
  sep_lines = sep_lines[c(-1, -length(reads))]
  
  tissue_name_full = NULL
  for(j in 1:(length(reads) - 1)){
    temp = rep(tissue_name[j], dim(reads[[j]][[2]])[1])
    tissue_name_full = c(tissue_name_full, temp)
  }
  
total_dim = 0
for(j in 1:(length(reads) - 1)){
  temp_dim = dim(reads[[j]][[2]])[1]
  total_dim = c(total_dim, total_dim[j] + temp_dim) 
}

res_pi_ordered = NULL
for(j in 1:(length(reads) - 1)){
  temp = res$pi[(total_dim[j]+1):total_dim[j+1], ]
  res_pi_ordered = rbind(res_pi_ordered, temp[order_reads[[j]], ])
}



  #     B = dim(res$phi)[2]
  #     J = floor(log2(B))
  #     
  #     phi_reflect = cbind(res$phi, res$phi[, B:(B - (2^(J + 1) - B - 1))])
  #     phi_reflect = cbind(phi_reflect, phi_reflect[, (2^(J+1):1)])
  #     
  
  res_phi_smooth = t(apply(res$phi, 1, ma, order = 10))
  y_max = max(res_phi_smooth, na.rm = TRUE)
  
  ind_pos = which(res$pi[, 3] > 0.9)
  ind_id = info_full[ind_pos]
  ind_col = factor(ind_id, labels = c(1, 5:10))
  ind_col = as.numeric(levels(ind_col))[ind_col]
  ind_col[ind_col == 7] = "orange"
  ind_col[ind_col == 8] = "gray47"
  ind_col[ind_col == 9] = "green4"
  ind_col[ind_col == 10] = "brown"


  pdf(paste0("cluster_res_all_", K, "_clusters_", gene_region[gi[g]], "_all.pdf"), height = 12, width = 14)
  par(mfrow = c(3, 1), mar = c(4, 4, 4, 1), oma = c(2, 2, 1, 0.5))
  plot(res_phi_smooth[1, ], ylim = c(0, y_max), type = 'l', col = 2, cex.lab = 1.2, xlab = "location (relative to TSS)", ylab = 'normalized intensity')
  for(k in 2:K){
    lines(res_phi_smooth[k, ], col = k + 1)
  }
  barplot(t(res$pi), col = 2:(K + 1), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
  axis(1, at = sep_lines_mid, labels = tissue_name, cex = 2, padj = -1, tick = FALSE)
  for(l in 1:length(ind_id)){
    axis(3, at = ind_pos[l], labels = ind_id[l], cex.axis = 0.7, col.axis = ind_col[l], las = 2)
  }
  mtext("tissue", 1, line = 2, cex = 1.2)
  mtext("membership proportion", 2, line = 4, cex = 1.2)
  abline(v = sep_lines, lwd = 2)

  plot(res_pi_ordered[, 2], ylim = c(0, 1), xaxt = 'n', xlab = "", ylab = "membership proportion", cex.lab = 1.5, cex.axis = 1.5, type = 'l', col = 3)
lines(res_pi_ordered[, 3], col = 4)
lines(res_pi_ordered[, 1], col = 2)
# lines(res_pi_ordered[, 4], col = 5)

axis(1, at = sep_lines_mid, labels = tissue_name, cex = 2, padj = -1, tick = FALSE)
mtext("tissue", 1, line = 2, cex = 1.2)
abline(v = sep_lines, lwd = 2)

  dev.off()


