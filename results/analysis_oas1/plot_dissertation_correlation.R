vcfdata = read.table("../../data/oas1/12.113000000-113800000.ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", header = FALSE)
vcfinfo = vcfdata[, 1:3]
snp_orig = vcfinfo[which(vcfinfo[, 3] == "rs10774671"), ]

genotype = read.table("../../data/oas1/oas1_yri_geno.012", header = FALSE)
indv = read.table("../../data/oas1/oas1_yri_geno.012.indv", header = FALSE)
pos = read.table("../../data/oas1/oas1_yri_geno.012.pos", header = FALSE)
snp_orig_pos = which(pos[, 2] == snp_orig[, 2])

row.names(genotype) = indv[, 1]

indv.id =  c("NA18486", "NA18498", "NA18499", "NA18501", "NA18504", "NA18505", "NA18507",
             "NA18508", "NA18510", "NA18511", "NA18516", "NA18517", "NA18519", "NA18520",
             "NA18522", "NA18523", "NA18852", "NA18855", "NA18858", "NA18861", "NA18862",
             "NA18870", "NA18871", "NA18909", "NA18912", "NA18913", "NA18916", "NA19093",
             "NA19098", "NA19099", "NA19101", "NA19102", "NA19108", "NA19114", "NA19116",
             "NA19119", "NA19127", "NA19128", "NA19130", "NA19131", "NA19137", "NA19138",
             "NA19140", "NA19143", "NA19144", "NA19147", "NA19152", "NA19153", "NA19159",
             "NA19160", "NA19171", "NA19172", "NA19190", "NA19192", "NA19193", "NA19200",
             "NA19201", "NA19203", "NA19204", "NA19209", "NA19210", "NA19222", "NA19225",
             "NA19238", "NA19239", "NA19257")

indv.id %in% row.names(genotype)

#Subset genotype matrix for 1000 Genomes individuals to those in the Pickrell paper
genotype = genotype[row.names(genotype) %in% indv.id, ]
genotype = genotype[, -1]

#####

load("../../results/analysis_oas1/res_smooth_3.RData")
res_smooth = res
load(paste0("../../data/oas1/OAS1.Robj"))

#Subset individuals in the Pickrell paper to those in 1000 Genomes(66 -> 55 individuals)
K = 3
g = Robj$g[indv.id %in% row.names(genotype)]
ordering = order(g)
smooth_pi_unordered = res_smooth$pi
smooth_pi_unordered_trunc = smooth_pi_unordered[indv.id %in% row.names(genotype), ]
smooth_pi_ordered_trunc = smooth_pi_unordered_trunc[ordering, ]

indv.id.ordered = indv.id[ordering]


#Construct #clusters by #snps correlation matrix, where (i,j)th entry is the correlation between cluster memberships (for the 55 samples) and genotypes (for the 55 samples)

cor.clust = matrix(0, nr = K, nc = dim(genotype)[2])
for(j in 1:dim(genotype)[2]){
  for(i in 1:K){
    cor.clust[i, j] = cor(smooth_pi_ordered_trunc[, i], genotype[ordering, j])
  }
}

col_cluster = c("green", "red", "blue")

pdf("association_smooth_3_cluster.pdf", height = 10, width = 8)
par(mfrow = c(K, 1), mar=c(5, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
for(i in 1:K){
  max_corr = which(abs(cor.clust[i, ]) == max(abs(cor.clust[i, ]), na.rm = TRUE))
  plot(cor.clust[i, ], pch = 16, cex = 0.7, xlab = "SNP index", ylab = "correlation", main = paste0("correlation between cluster memberships and genotypes for the ", col_cluster[i], " cluster"))
  points(max_corr, cor.clust[i, max_corr], pch = 16, cex = 1.5, col = 4)
  abline(v = snp_orig_pos, lty = 2, col = 2)
}
dev.off()




load("../../results/analysis_oas1/res_nosmooth_3.RData")
res_nosmooth = res
load(paste0("../../data/oas1/OAS1.Robj"))

#Subset individuals in the Pickrell paper to those in 1000 Genomes(66 -> 55 individuals)
K = 3
g = Robj$g[indv.id %in% row.names(genotype)]
ordering = order(g)
nosmooth_pi_unordered = res_nosmooth$pi[, c(2, 1, 3)]
nosmooth_pi_unordered_trunc = nosmooth_pi_unordered[indv.id %in% row.names(genotype), ]
nosmooth_pi_ordered_trunc = nosmooth_pi_unordered_trunc[ordering, ]

indv.id.ordered = indv.id[ordering]


#Construct #clusters by #snps correlation matrix, where (i,j)th entry is the correlation between cluster memberships (for the 55 samples) and genotypes (for the 55 samples)

cor.clust = matrix(0, nr = K, nc = dim(genotype)[2])
for(j in 1:dim(genotype)[2]){
  for(i in 1:K){
    cor.clust[i, j] = cor(nosmooth_pi_ordered_trunc[, i], genotype[ordering, j])
  }
}

col_cluster = c("green", "red", "blue")

pdf("association_unsmooth_3_cluster.pdf", height = 10, width = 8)
par(mfrow = c(K, 1), mar=c(5, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
for(i in 1:K){
  max_corr = which(abs(cor.clust[i, ]) == max(abs(cor.clust[i, ]), na.rm = TRUE))
  plot(cor.clust[i, ], pch = 16, cex = 0.7, xlab = "SNP index", ylab = "correlation", main = paste0("correlation between cluster memberships and genotypes for the ", col_cluster[i], " cluster"))
  points(max_corr, cor.clust[i, max_corr], pch = 16, cex = 1.5, col = 4)
  abline(v = snp_orig_pos, lty = 2, col = 2)
}
dev.off()





load("../../results/analysis_oas1/res_smooth_2.RData")
res_smooth = res
load(paste0("../../data/oas1/OAS1.Robj"))

#Subset individuals in the Pickrell paper to those in 1000 Genomes(66 -> 55 individuals)
K = 2
g = Robj$g[indv.id %in% row.names(genotype)]
ordering = order(g)
smooth_pi_unordered = res_smooth$pi[, c(2, 1)]
smooth_pi_unordered_trunc = smooth_pi_unordered[indv.id %in% row.names(genotype), ]
smooth_pi_ordered_trunc = smooth_pi_unordered_trunc[ordering, ]

indv.id.ordered = indv.id[ordering]


#Construct #clusters by #snps correlation matrix, where (i,j)th entry is the correlation between cluster memberships (for the 55 samples) and genotypes (for the 55 samples)

cor.clust = matrix(0, nr = K, nc = dim(genotype)[2])
for(j in 1:dim(genotype)[2]){
  for(i in 1:K){
    cor.clust[i, j] = cor(smooth_pi_ordered_trunc[, i], genotype[ordering, j])
  }
}

col_cluster = c("green", "red")

pdf("association_smooth_2_cluster.pdf", height = 7, width = 8)
par(mfrow = c(K, 1), mar=c(5, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
for(i in 1:K){
  max_corr = which(abs(cor.clust[i, ]) == max(abs(cor.clust[i, ]), na.rm = TRUE))
  plot(cor.clust[i, ], pch = 16, cex = 0.7, xlab = "SNP index", ylab = "correlation", main = paste0("correlation between cluster memberships and genotypes for the ", col_cluster[i], " cluster"))
  points(max_corr, cor.clust[i, max_corr], pch = 16, cex = 1.5, col = 4)
  abline(v = snp_orig_pos, lty = 2, col = 2)
}
dev.off()




load("../../results/analysis_oas1/res_nosmooth_2.RData")
res_nosmooth = res
load(paste0("../../data/oas1/OAS1.Robj"))

#Subset individuals in the Pickrell paper to those in 1000 Genomes(66 -> 55 individuals)
K = 2
g = Robj$g[indv.id %in% row.names(genotype)]
ordering = order(g)
nosmooth_pi_unordered = res_nosmooth$pi[, c(2, 1)]
nosmooth_pi_unordered_trunc = nosmooth_pi_unordered[indv.id %in% row.names(genotype), ]
nosmooth_pi_ordered_trunc = nosmooth_pi_unordered_trunc[ordering, ]

indv.id.ordered = indv.id[ordering]


#Construct #clusters by #snps correlation matrix, where (i,j)th entry is the correlation between cluster memberships (for the 55 samples) and genotypes (for the 55 samples)

cor.clust = matrix(0, nr = K, nc = dim(genotype)[2])
for(j in 1:dim(genotype)[2]){
  for(i in 1:K){
    cor.clust[i, j] = cor(nosmooth_pi_ordered_trunc[, i], genotype[ordering, j])
  }
}

col_cluster = c("green", "red")

pdf("association_unsmooth_2_cluster.pdf", height = 7, width = 8)
par(mfrow = c(K, 1), mar=c(5, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
for(i in 1:K){
  max_corr = which(abs(cor.clust[i, ]) == max(abs(cor.clust[i, ]), na.rm = TRUE))
  plot(cor.clust[i, ], pch = 16, cex = 0.7, xlab = "SNP index", ylab = "correlation", main = paste0("correlation between cluster memberships and genotypes for the ", col_cluster[i], " cluster"))
  points(max_corr, cor.clust[i, max_corr], pch = 16, cex = 1.5, col = 4)
  abline(v = snp_orig_pos, lty = 2, col = 2)
}
dev.off()






load("../../results/analysis_oas1/res_smooth_4.RData")
res_smooth = res
load(paste0("../../data/oas1/OAS1.Robj"))

#Subset individuals in the Pickrell paper to those in 1000 Genomes(66 -> 55 individuals)
K = 4
g = Robj$g[indv.id %in% row.names(genotype)]
ordering = order(g)
smooth_pi_unordered = res_smooth$pi
smooth_pi_unordered_trunc = smooth_pi_unordered[indv.id %in% row.names(genotype), ]
smooth_pi_ordered_trunc = smooth_pi_unordered_trunc[ordering, ]

indv.id.ordered = indv.id[ordering]


#Construct #clusters by #snps correlation matrix, where (i,j)th entry is the correlation between cluster memberships (for the 55 samples) and genotypes (for the 55 samples)

cor.clust = matrix(0, nr = K, nc = dim(genotype)[2])
for(j in 1:dim(genotype)[2]){
  for(i in 1:K){
    cor.clust[i, j] = cor(smooth_pi_ordered_trunc[, i], genotype[ordering, j])
  }
}

col_cluster = c("green", "red", "blue", "cyan")

pdf("association_smooth_4_cluster.pdf", height = 13, width = 8)
par(mfrow = c(K, 1), mar=c(5, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
for(i in 1:K){
  max_corr = which(abs(cor.clust[i, ]) == max(abs(cor.clust[i, ]), na.rm = TRUE))
  plot(cor.clust[i, ], pch = 16, cex = 0.7, xlab = "SNP index", ylab = "correlation", main = paste0("correlation between cluster memberships and genotypes for the ", col_cluster[i], " cluster"))
  points(max_corr, cor.clust[i, max_corr], pch = 16, cex = 1.5, col = 4)
  abline(v = snp_orig_pos, lty = 2, col = 2)
}
dev.off()






load("../../results/analysis_oas1/res_nosmooth_4.RData")
res_nosmooth = res
load(paste0("../../data/oas1/OAS1.Robj"))

#Subset individuals in the Pickrell paper to those in 1000 Genomes(66 -> 55 individuals)
K = 4
g = Robj$g[indv.id %in% row.names(genotype)]
ordering = order(g)
nosmooth_pi_unordered = res_nosmooth$pi[, c(4, 1, 3, 2)]
nosmooth_pi_unordered_trunc = nosmooth_pi_unordered[indv.id %in% row.names(genotype), ]
nosmooth_pi_ordered_trunc = nosmooth_pi_unordered_trunc[ordering, ]

indv.id.ordered = indv.id[ordering]


#Construct #clusters by #snps correlation matrix, where (i,j)th entry is the correlation between cluster memberships (for the 55 samples) and genotypes (for the 55 samples)

cor.clust = matrix(0, nr = K, nc = dim(genotype)[2])
for(j in 1:dim(genotype)[2]){
  for(i in 1:K){
    cor.clust[i, j] = cor(nosmooth_pi_ordered_trunc[, i], genotype[ordering, j])
  }
}

col_cluster = c("green", "red", "blue", "cyan")

pdf("association_unsmooth_4_cluster.pdf", height = 13, width = 8)
par(mfrow = c(K, 1), mar=c(5, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
for(i in 1:K){
  max_corr = which(abs(cor.clust[i, ]) == max(abs(cor.clust[i, ]), na.rm = TRUE))
  plot(cor.clust[i, ], pch = 16, cex = 0.7, xlab = "SNP index", ylab = "correlation", main = paste0("correlation between cluster memberships and genotypes for the ", col_cluster[i], " cluster"))
  points(max_corr, cor.clust[i, max_corr], pch = 16, cex = 1.5, col = 4)
  abline(v = snp_orig_pos, lty = 2, col = 2)
}
dev.off()
