library(forecast)

gene_region = read.table("../../supplemental/gtex/gene_region.txt", stringsAsFactors = FALSE)
gene_region = unlist(gene_region, use.names = FALSE)
gene_region = gsub("\\:", "_", gene_region)
gene_region = gsub("\\-", "_", gene_region)


for(K in 2:5){
  for(g in 1:length(gene_region)){
    load_name = paste("cluster_res_all", K, "clusters", gene_region[g], sep = "_")
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
    
#     B = dim(res$phi)[2]
#     J = floor(log2(B))
#     
#     phi_reflect = cbind(res$phi, res$phi[, B:(B - (2^(J + 1) - B - 1))])
#     phi_reflect = cbind(phi_reflect, phi_reflect[, (2^(J+1):1)])
#     

    res_phi_smooth = t(apply(res$phi, 1, ma, order = 10))
    y_max = max(res_phi_smooth, na.rm = TRUE)

    pdf(paste0("cluster_res_all_", K, "_clusters_", gene_region[g], ".pdf"), height = 8, width = 12)
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(2, 2, 0.5, 0.5))
    plot(res_phi_smooth[1, ], ylim = c(0, y_max), type = 'l', col = 2, cex.lab = 1.2, xlab = "location (relative to TSS)", ylab = 'normalized intensity')
    for(k in 2:K){
      lines(res_phi_smooth[k, ], col = k + 1)
    }
    barplot(t(res$pi), col = 2:(K + 1), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
    axis(1, at = sep_lines_mid, labels = tissue_name, cex = 2, padj = -1, tick = FALSE)
    mtext("tissue", 1, line = 2, cex = 1.2)
    mtext("membership proportion", 2, line = 4, cex = 1.2)
    abline(v = sep_lines, lwd = 2)
    dev.off()
  }
}
