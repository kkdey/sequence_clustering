dir.name = "~/projects/sequence_clustering"
library(multiseq)
lreg = 2^10
peaks = read.table(file.path(dir.name, "data", "k562ctcf", "wgEncodeBroadHistoneK562CtcfStdAlnRep0.bed"))
samplesheet = file.path(dir.name, "src/analysis_k562ctcf", "samplesheet_k562ctcf")
peak.chr = peaks[, 1]
peak.center = ceiling((peaks[, 2] + peaks[, 3])/2)
peak.start = peak.center - lreg/2
peak.end = peak.start + lreg - 1
nreg = length(peak.chr)

data.sig = matrix(0, nreg, lreg)
for(i in 1:nreg){
	region = paste0(peak.chr[i], ":", peak.start[i], "-", peak.end[i])
	data.sig[i, ] = get.counts(samplesheet, region)[1,]
	print(i)
}
save(data.sig, peak.chr, peak.center, peak.start, peak.start, file = file.path(dir.name, "data", "k562ctcf", "data.sig.Robj"))
