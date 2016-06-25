library(gplots)

load("res_nosmooth_3.RData")
res_nosmooth = res
load("res_smooth_3.RData")
res_smooth = res

pdf("cluster_mean_3_cluster.pdf", width = 6, height = 4)
par(mfrow = c(2, 1), mar = c(1, 4, 1, 1), oma = c(2, 2, 0.5, 0.5))
plot(res_nosmooth$phi[1, ], xlab = "", ylab = "normalized intensity", xaxt = "n", col = 2, ylim = c(0, 0.02), type = 'l')
lines(res_nosmooth$phi[2, ], col = 3)
lines(res_nosmooth$phi[3, ], col = 4)
plot(res_smooth$phi[1, ], xlab = "", ylab = "normalized intensity", xaxt = "n", col = 3, ylim = c(0, 0.02), type = 'l')
lines(res_smooth$phi[2, ], col = 2)
lines(res_smooth$phi[3, ], col = 4)
dev.off()

K = dim(res_nosmooth$pi)[2]
nosmooth_pi_unordered = res_nosmooth$pi
ordering = order(Robj$g)
nosmooth_pi_ordered = nosmooth_pi_unordered[ordering, c(2, 1, 3)]
smooth_pi_unordered = res_smooth$pi
smooth_pi_ordered = smooth_pi_unordered[ordering, ]
sep_lines = cumsum(table(Robj$g))[1:2]
sep_lines_ext = c(0, sep_lines, dim(res_nosmooth$pi)[1])
sep_lines_mid = (sep_lines_ext[1:3] + sep_lines_ext[2:4])/2

pdf("admixture_3_cluster.pdf", width = 12, height = 8)
par(mfrow = c(2, 1), mar=c(2, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
barplot(t(nosmooth_pi_ordered), col = c(3, 2, 4:(K + 1)), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
axis(1, at = sep_lines_mid, labels = 0:2, cex = 2, padj = -1, tick = FALSE)
mtext("genotype", 1, line = 2, cex = 2)
mtext("membership proportion", 2, line = 4, cex = 1.2)
abline(v = sep_lines, lwd = 2)
barplot(t(smooth_pi_ordered), col = c(3, 2, 4:(K + 1)), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
axis(1, at = sep_lines_mid, labels = 0:2, cex = 2, padj = -1, tick = FALSE)
mtext("genotype", 1, line = 2, cex = 2)
mtext("membership proportion", 2, line = 4, cex = 1.2)
abline(v = sep_lines, lwd = 2)
dev.off()



col = c(rgb(seq(0, 1, length = 15), 1, seq(0, 1, length = 15)), rgb(1, seq(1, 0 ,length = 15), seq(1, 0, length = 15)))

pdf("heatmap_unsmooth_3_cluster.pdf", height = 8, width = 12)
heatmap.2(cor(t(res_nosmooth$pi), method = "pearson"), cexRow = 0.8, cexCol = 1.15, labCol = Robj$g, labRow = Robj$g, scale = "none", trace = "none", distfun = function(x) dist(x, method = "euclidean"), col = col, hclustfun = function(x) hclust(x,method="average"), breaks = c(seq(-1, 0, 0.1), seq(0.05, 1, length.out = 20)))
dev.off()

pdf("heatmap_smooth_3_cluster.pdf", height = 8, width = 12)
heatmap.2(cor(t(res_smooth$pi), method = "pearson"), cexRow = 0.8, cexCol = 1.15, labCol = Robj$g, labRow = Robj$g, scale = "none", trace = "none", distfun = function(x) dist(x, method = "euclidean"), col = col, hclustfun = function(x) hclust(x,method="average"), breaks = c(seq(-1, 0, 0.1), seq(0.05, 1, length.out = 20)))
dev.off()



#############



load("res_nosmooth_2.RData")
res_nosmooth = res
load("res_smooth_2.RData")
res_smooth = res

pdf("cluster_mean_2_cluster.pdf", width = 6, height = 4)
par(mfrow = c(2, 1), mar = c(1, 4, 1, 1), oma = c(2, 2, 0.5, 0.5))
plot(res_nosmooth$phi[1, ], xlab = "", ylab = "normalized intensity", xaxt = "n", col = 2, ylim = c(0, 0.02), type = 'l')
lines(res_nosmooth$phi[2, ], col = 3)
plot(res_smooth$phi[1, ], xlab = "", ylab = "normalized intensity", xaxt = "n", col = 2, ylim = c(0, 0.02), type = 'l')
lines(res_smooth$phi[2, ], col = 3)
dev.off()

K = dim(res_nosmooth$pi)[2]
nosmooth_pi_unordered = res_nosmooth$pi[, c(2, 1)]
ordering = order(Robj$g)
nosmooth_pi_ordered = nosmooth_pi_unordered[ordering, ]
smooth_pi_unordered = res_smooth$pi[, c(2, 1)]
smooth_pi_ordered = smooth_pi_unordered[ordering, ]
sep_lines = cumsum(table(Robj$g))[1:2]
sep_lines_ext = c(0, sep_lines, dim(res_nosmooth$pi)[1])
sep_lines_mid = (sep_lines_ext[1:3] + sep_lines_ext[2:4])/2

pdf("admixture_2_cluster.pdf", width = 12, height = 8)
par(mfrow = c(2, 1), mar=c(2, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
barplot(t(nosmooth_pi_ordered), col = c(3, 2), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
axis(1, at = sep_lines_mid, labels = 0:2, cex = 2, padj = -1, tick = FALSE)
mtext("genotype", 1, line = 2, cex = 2)
mtext("membership proportion", 2, line = 4, cex = 1.2)
abline(v = sep_lines, lwd = 2)
barplot(t(smooth_pi_ordered), col = c(3, 2), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
axis(1, at = sep_lines_mid, labels = 0:2, cex = 2, padj = -1, tick = FALSE)
mtext("genotype", 1, line = 2, cex = 2)
mtext("membership proportion", 2, line = 4, cex = 1.2)
abline(v = sep_lines, lwd = 2)
dev.off()



col = c(rgb(seq(0, 1, length = 15), 1, seq(0, 1, length = 15)), rgb(1, seq(1, 0 ,length = 15), seq(1, 0, length = 15)))

pdf("heatmap_unsmooth_2_cluster.pdf", height = 8, width = 12)
heatmap.2(cor(t(res_nosmooth$pi), method = "pearson"), cexRow = 0.8, cexCol = 1.15, labCol = Robj$g, labRow = Robj$g, scale = "none", trace = "none", distfun = function(x) dist(x, method = "euclidean"), col = col, hclustfun = function(x) hclust(x,method="average"), breaks = c(seq(-1, 0, 0.1), seq(0.05, 1, length.out = 20)))
dev.off()

pdf("heatmap_smooth_2_cluster.pdf", height = 8, width = 12)
heatmap.2(cor(t(res_smooth$pi), method = "pearson"), cexRow = 0.8, cexCol = 1.15, labCol = Robj$g, labRow = Robj$g, scale = "none", trace = "none", distfun = function(x) dist(x, method = "euclidean"), col = col, hclustfun = function(x) hclust(x,method="average"), breaks = c(seq(-1, 0, 0.1), seq(0.05, 1, length.out = 20)))
dev.off()


##############

load("res_nosmooth_4.RData")
res_nosmooth = res
load("res_smooth_4.RData")
res_smooth = res

pdf("cluster_mean_4_cluster.pdf", width = 6, height = 4)
par(mfrow = c(2, 1), mar = c(1, 4, 1, 1), oma = c(2, 2, 0.5, 0.5))
plot(res_nosmooth$phi[1, ], xlab = "", ylab = "normalized intensity", xaxt = "n", col = 2, ylim = c(0, 0.02), type = 'l')
lines(res_nosmooth$phi[2, ], col = 5)
lines(res_nosmooth$phi[3, ], col = 4)
lines(res_nosmooth$phi[4, ], col = 3)
plot(res_smooth$phi[1, ], xlab = "", ylab = "normalized intensity", xaxt = "n", col = 3, ylim = c(0, 0.02), type = 'l')
lines(res_smooth$phi[2, ], col = 2)
lines(res_smooth$phi[3, ], col = 4)
lines(res_smooth$phi[4, ], col = 5)
dev.off()

K = dim(res_nosmooth$pi)[2]
nosmooth_pi_unordered = res_nosmooth$pi[, c(4, 1, 3, 2)]
ordering = order(Robj$g)
nosmooth_pi_ordered = nosmooth_pi_unordered[ordering, ]
smooth_pi_unordered = res_smooth$pi
smooth_pi_ordered = smooth_pi_unordered[ordering, ]
sep_lines = cumsum(table(Robj$g))[1:2]
sep_lines_ext = c(0, sep_lines, dim(res_nosmooth$pi)[1])
sep_lines_mid = (sep_lines_ext[1:3] + sep_lines_ext[2:4])/2

pdf("admixture_4_cluster.pdf", width = 12, height = 8)
par(mfrow = c(2, 1), mar=c(2, 4, 2, 2), oma = c(2, 2, 0.2, 0.2))
barplot(t(nosmooth_pi_ordered), col = c(3, 2, 4:(K + 1)), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
axis(1, at = sep_lines_mid, labels = 0:2, cex = 2, padj = -1, tick = FALSE)
mtext("genotype", 1, line = 2, cex = 2)
mtext("membership proportion", 2, line = 4, cex = 1.2)
abline(v = sep_lines, lwd = 2)
barplot(t(smooth_pi_ordered), col = c(3, 2, 4:(K + 1)), axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)
axis(1, at = sep_lines_mid, labels = 0:2, cex = 2, padj = -1, tick = FALSE)
mtext("genotype", 1, line = 2, cex = 2)
mtext("membership proportion", 2, line = 4, cex = 1.2)
abline(v = sep_lines, lwd = 2)
dev.off()



col = c(rgb(seq(0, 1, length = 15), 1, seq(0, 1, length = 15)), rgb(1, seq(1, 0 ,length = 15), seq(1, 0, length = 15)))

pdf("heatmap_unsmooth_4_cluster.pdf", height = 8, width = 12)
heatmap.2(cor(t(res_nosmooth$pi), method = "pearson"), cexRow = 0.8, cexCol = 1.15, labCol = Robj$g, labRow = Robj$g, scale = "none", trace = "none", distfun = function(x) dist(x, method = "euclidean"), col = col, hclustfun = function(x) hclust(x,method="average"), breaks = c(seq(-1, 0, 0.1), seq(0.05, 1, length.out = 20)))
dev.off()

pdf("heatmap_smooth_4_cluster.pdf", height = 8, width = 12)
heatmap.2(cor(t(res_smooth$pi), method = "pearson"), cexRow = 0.8, cexCol = 1.15, labCol = Robj$g, labRow = Robj$g, scale = "none", trace = "none", distfun = function(x) dist(x, method = "euclidean"), col = col, hclustfun = function(x) hclust(x,method="average"), breaks = c(seq(-1, 0, 0.1), seq(0.05, 1, length.out = 20)))
dev.off()

