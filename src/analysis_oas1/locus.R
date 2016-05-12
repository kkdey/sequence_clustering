#usage
Rscript locus.R YUR chr12 113344739 113369990 0 2 rs10774671 112

#chr="chr12"
#data="YUR"
#locus.start=113344739
#locus.end=113369990
#snp="rs10774671"
#list.loci.line=112
#bin.size=2
#min.reads=1

course.repodir    <- scan(".course.repodir.txt", what=character())
ash.repodir       <- scan(".ash.repodir.txt", what=character())
RNAsplice.repodir <- scan(".RNAsplice.repodir.txt", what= character())
source(file.path(RNAsplice.repodir, "code/R/PoissonBinomialpairedend.funcs.R"))
source(file.path(RNAsplice.repodir, "code/R/utils.R"))
source(file.path(course.repodir, "/code/Rcode/get_shrunken_WCs.R"))
setwd(RNAsplice.repodir)



#********************************************************************
#
#     get arguments
#
#********************************************************************
args        <- commandArgs(TRUE)
data        <- args[1] 
chr         <- args[2]
locus.start <- as.numeric(args[3])
locus.end   <- as.numeric(args[4])
min.reads   <- as.numeric(args[5])
bin.size    <- as.numeric(args[6])
force.locus <- TRUE
#min.percent <- 0



#**************************************************************************
#
#      define data
#
#**************************************************************************
if (data=="YUR"){
    source.file    <- file.path(RNAsplice.repodir, "/code/scripts/utilities.sh")
    snp            <- args[7]
    list.loci.line <- args[8]
    plotListLociLine <- function(locus.start, locus.end, file, line){
        DiffLoci   <- data.frame(lapply(read.table(file, fill=1, header=0 ), as.character), stringsAsFactors=FALSE)
        sqtl.start <- as.numeric(DiffLoci[line,4])
        sqtl.end   <- as.numeric(DiffLoci[line,5])
        rect(sqtl.start, 0, sqtl.end, 1, col=rgb(0,1,0,0.5), border=NA)
    }
    #edited
    individuals <- c("NA18486","NA18498","NA18499","NA18501","NA18502","NA18504_2","NA18505_2","NA18507","NA18508","NA18510","NA18511","NA18516_2","NA18517","NA18519","NA18520","NA18522","NA18523","NA18852","NA18855","NA18856_2","NA18858","NA18861","NA18862","NA18870","NA18871","NA18909","NA18912_2","NA18913","NA18916","NA19093","NA19098","NA19099","NA19101","NA19102","NA19108","NA19114","NA19116","NA19119","NA19127","NA19128","NA19130","NA19131","NA19137","NA19138","NA19140","NA19143","NA19144","NA19147","NA19152","NA19153","NA19159","NA19160","NA19171_2","NA19172","NA19190","NA19192","NA19193","NA19200","NA19201","NA19203","NA19204","NA19209","NA19210","NA19222","NA19225","NA19238","NA19239","NA19257")
    individuals.id <- c("NA18486","NA18498","NA18499","NA18501","NA18502","NA18504","NA18505","NA18507","NA18508","NA18510","NA18511","NA18516","NA18517","NA18519","NA18520","NA18522","NA18523","NA18852","NA18855","NA18856","NA18858","NA18861","NA18862","NA18870","NA18871","NA18909","NA18912","NA18913","NA18916","NA19093","NA19098","NA19099","NA19101","NA19102","NA19108","NA19114","NA19116","NA19119","NA19127","NA19128","NA19130","NA19131","NA19137","NA19138","NA19140","NA19143","NA19144","NA19147","NA19152","NA19153","NA19159","NA19160","NA19171","NA19172","NA19190","NA19192","NA19193","NA19200","NA19201","NA19203","NA19204","NA19209","NA19210","NA19222","NA19225","NA19238","NA19239","NA19257")
    bamfiles           <- paste0("/mnt/lustre/home/epantaleo/data/YUR/YUR_tophat/YUR_tophat-2.0.0/", individuals, "/accepted_hits.bam")
    string.individuals <- paste(individuals, collapse=" ")
    pcs                <- read.table("/mnt/lustre/home/surbut/src/RNAsplice/analysis/single_end/pc_loadings")
    pcs                <- as.matrix(pcs[individuals.id, 1:10])
    command            <- paste("source", source.file, "; individuals=(", string.individuals, "); write_genotypes_YUR individuals[@]", snp, chr)
    genotype           <- system(command, intern = TRUE)
    genotype           <- genotype[1:(length(genotype)-1)]

    #genotype[genotype=="NN"] <- NA
    index              <- which(genotype!="NN")
    genotype           <- genotype[index]
    g                  <- factor(genotype, exclude = NA)
    g                  <- match(g, levels(g)) - 1

    qc.file            <- "/mnt/lustre/home/epantaleo/data/YUR/YUR_tophat/YUR_tophat-2.0.0/qc_YUR_tophat-2.0.0.gz"
    command            <- paste("source", source.file, "; individuals=(", string.individuals, "); get_totalcount_YUR ", qc.file, "individuals[@]")
    read.depth         <- as.numeric(unlist(strsplit(system(command, intern = TRUE),",")))
    read.depth         <- read.depth[index]
    bamfiles           <- bamfiles[index]
    pcs                <- pcs[index,]
    individuals        <- individuals[index]
    force.locus        <- TRUE
#    merged.bed         <- 
    list.loci.file     <- "/mnt/lustre/home/epantaleo/scripts/list_loci_hg19_all"
}

#*************************************************************************
#
#     get M.0 (L by n.col)
#     L is number of individuals (68 in the case of YUR data)
#     n.col=locus.end - locus.start + 1
#     elements of M.0 are number of reads that start at each position in [locus.start,locus.end]
#     in each individual
#
#************************************************************************ 
#extend locus
if (force.locus==FALSE){
    extended.locus <- extend.locus(chr, locus.start, locus.end, merged.bed)
    locus.start    <- extended.locus$start
    locus.end      <- extended.locus$end
    print(paste("Locus has been extended from [",
                locus.start,
                ",",
                locus.end,
                "] to [",
                extended.locus$start,
                ",",
                extended.locus$end,
                "]"))
    locus.start    <- extended.locus$start
    locus.end      <- extended.locus$end
}
#M.0 is a sparse matrix (see R package "Matrix")
print("get data")
system.time(M.0  <- get.1D.Matrix(bamfiles, chr, locus.start, locus.end))
N.0  <- ncol(M.0)
L    <- nrow(M.0)

#*************************************************************************
#
#     filter M.0 to get M (L by N)
#     where N is a power of 2
#
#************************************************************************  
#note: function reduce.size.Matrix creates a matrix that has number of col a power of 2
#(after filtering)
#it returns error if number of columns after filtering is greater than 2^20
print("change size of data")
M          <- reduce.size.Matrix(M.0, min.reads, bin.size)
subset.col <- M$cols.filter
N          <- ncol(M$M)
J          <- log2(N)
T          <- log2(N)*N

locus.name <- paste(chr, locus.start, locus.end, sep=".")
dir.name   <- file.path(RNAsplice.repodir, "analysis", "single_end", "results")
dir.create(dir.name)
dir.name   <- file.path(dir.name, locus.name)
dir.create(dir.name)
dir.name   <- file.path(dir.name, paste(min.reads, bin.size, sep="_"))
dir.create(dir.name)

###Convert the shrunken wavelet coefficients to a matrix suitable for PCA###
print("get shrunk wavelet coefficients in a TI table")
ptm        <- proc.time()
alpha.mean <- matrix(0, nrow=T, ncol=L)
mu.mean    <- NULL
good.people <-NULL
for (i in 1:L){
    m              <- M$M[i,]
    shrunkenWCs    <- shrunken_WCs(m, read.depth[i])
    if(!is.null(shrunkenWCs)){
    	alpha.mean[,i] <- as.vector(shrunkenWCs$alpha.mean)
    	mu.mean        <- c(mu.mean, shrunkenWCs$mu.mean)
    	   good.guys=i
    	   good.people=c(good.people,good.guys)
    	   }
   }
good.people=matrix(good.people)
write.table(good.people, file=file.path(dir.name, "good_people.txt"), row.names=FALSE, col.names=FALSE, append = FALSE)

alpha.Mean <- rowMeans(alpha.mean)
print(proc.time() - ptm)
save(list=ls(),file=file.path(dir.name, paste0(locus.name,".Robj")))

print("perform svd")
                                        #arbitrary choice of number of components K
if (L < 15) K <- 2 else K <- 10
pca <- svd(alpha.mean - alpha.Mean)
D   <- pca$d
U   <- pca$u
V   <- pca$v

###K=which(cumsum(D^2/sum(D^2))>=0.85)[1]

###show that this makes a big difference as opposed to performing PCA on simply the data matrix##
pca2 <- svd(t(M$M - colMeans(M$M)))
D2   <- pca2$d
U2   <- pca2$u
V2   <- pca2$v

#cumsum(D^2/sum(D^2))

#for(K in 1:L){
#manova.full.2  <- lm(g ~ V[,1:K] + mu.mean + pcs)
#manova.small.2 <- lm(g ~ mu.mean + pcs)
#ppval2.b[K]         <- anova(manova.full.2, manova.small.2)[[2,6]] 
#}

#ppval2.b=na.omit(ppval2.b)
#K= which(ppval2.b==min(ppval2.b))

pmat      <- matrix(NA,nrow=6,ncol=1) 

mu.full   <- lm(mu.mean ~ g + pcs) 
ppval1    <- summary(mu.full)$coefficients[2,4] 
pmat[1,]  <- ppval1

manova.full.2  <- lm(g ~ V[,1:K] + mu.mean + pcs)
manova.small.2 <- lm(g ~ mu.mean + pcs)
ppval2         <- anova(manova.full.2, manova.small.2)[[2,6]] 
pmat[2,]       <- ppval2

manova.full.3  <- lm(g ~ V[,1:K])
manova.small.3 <- lm(g ~ 1)
ppval3         <- anova(manova.full.3, manova.small.3)[[2,6]]
pmat[3,]       <- ppval3

manova.full.4  <- lm(g ~ V[,1:K] + mu.mean + pcs)
manova.small.4 <- lm(g ~ pcs)
ppval4         <- anova(manova.full.4, manova.small.4)[[2,6]]
pmat[4,]       <- ppval4

manova.full.5  <- lm(g ~ V[,1:K] + mu.mean)
manova.small.5 <- lm(g ~ 1)
ppval5         <- anova(manova.full.5, manova.small.5)[[2,6]]
pmat[5,]       <- ppval5 


manova.full.6  <- lm(g ~ V[,1:K] + mu.mean+pcs)
manova.small.6 <- lm(g ~ pcs)
ppval6         <- anova(manova.full.6, manova.small.6)[[2,6]]
pmat[6,]       <- ppval6 

write.table(pmat, file=file.path(dir.name, "pval.txt"), row.names=FALSE, col.names=FALSE, append = FALSE)


                                        #****************************************************************
                                        #
                                        #     make some plots
                                        #
                                        #****************************************************************
    
#pdf(file=file.path(dir.name, "ProportionVariationExplained.pdf"))
#par(mfrow=c(2,1))
#plot(D^2/sum(D^2), main="ProportionVarianceExplained with WC", xlab="PC")
#plot(D2^2/sum(D2^2), main="ProportionVarianceExplained w/o WC", xlab="PC")
#dev.off()

test1  <- lm(g~V[,1:K] + mu.mean) 
summary(test1) 
manova <- lm(g~V[,1:K] + mu.mean)
summary(manova) 
#test1 <- lm(g~V2[,1:K])
#summary(test1)

save(list=ls(),file=file.path(dir.name, paste0(locus.name,".Robj")))

pdf(file=file.path(dir.name, "IndividualsbyEigenWC.pdf"))
par(mfrow=c(2,1))
plot(V[,1], V[,2], main="with WC", pch=g+2, col=g+2)
plot(V2[,1], V2[,2], main="w/o WC", pch=g+2, col=g+2)
dev.off()

###Look at eigencurves in data space, first look at original curves##
##Transform back into the original data space and then add back in eliminated rows###
eigencurves         <- matrix(NA, nrow=N, ncol=K)
eigendeviations     <- matrix(NA, nrow=N, ncol=K)
for (k in 1:K){
    #eigencurves[,k]    <- reverseWT.PCA(matrix(U[,k]+alpha.Mean, nrow=J, byrow=FALSE), 0.025)$est.sm
    eigencurves[,k]     <- reverseWT.PCA(matrix(U[,k] + alpha.Mean, nrow=J, byrow=FALSE), threshold=NULL, highlevelfree=2, logscale=FALSE)$est
    eigendeviations[,k] <- reverseWT.PCA(matrix(U[,k], nrow=J, byrow=FALSE), threshold=NULL, highlevelfree=2, logscale=FALSE)$est
}
save(list=ls(),file=file.path(dir.name, paste0(locus.name,".Robj")))

pdf(file=file.path(dir.name, "eigendeviation_curves.in.data.space.pdf"))
ymin <- min(eigendeviations)
ymax <- max(eigendeviations)
par(mfrow=c(5,1))
plot(as.matrix(eigendeviations[,1]), ylab="y", type="l", ylim=c(ymin,ymax), col="blue", main="eigendeviation 1 (blue), 2 (red), 3 (green)")
points(as.matrix(eigendeviations[,2]), ylab="y", type="l", col="red")
points(as.matrix(eigendeviations[,3]), ylab="y", type="l", col="green")
plot(as.vector(M$M[1,]), type="l", ylab="# reads", xlab=NULL,  main="original data: first individual")
plot(reconstruct.1D.Matrix(t(as.matrix(eigendeviations[,1])), N.0, subset.col, bin.size), ylab="y", pch=".", ylim=c(ymin,ymax),  xlab=NULL, col="blue", main="eigendeviation 1 in original space")
points(reconstruct.1D.Matrix(t(as.matrix(eigendeviations[,2])), N.0, subset.col, bin.size), pch=".", col="red")
points(reconstruct.1D.Matrix(t(as.matrix(eigendeviations[,3])), N.0, subset.col, bin.size), pch=".", col="green")
plot(reconstruct.1D.Matrix(t(as.matrix(M$M[1,])), N.0, subset.col, bin.size), ylab="# reads", type="l",  xlab=NULL, main="first individual in original space (original data)")

hg19.annotation.file.gz <- "/mnt/lustre/home/epantaleo/data/annotations/hg19.ensGene.gp.gz"
Transcripts             <- get.Transcripts(hg19.annotation.file.gz, chr, locus.start, locus.end)
plotGenePred(Transcripts, locus.start, locus.end)
if (data=="YUR")
    plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)
dev.off()

                                        #******************************************************************
                                        #
                                        #     get p
                                        #
                                        #******************************************************************
                                        #p_5 is the probability of using a 5' splice site (ss) when encountered
                                        #p_3 is the probability of jumping to a 3'ss given that it is considered
                                        #psi=logit(p)
                                        #p_5=s_5/(s_5+c_5)
                                        #s_5=number of reads that use the specific 5'ss
                                        #c_5=number of reads that map to the base following the 5'ss

print("compute mvdat,tcut=0,no.total.intensity")
tcut       <- 0
pval   <- rep(NA, K+1)
mvdat  <- matrix(NA, nr = 4, nc = K)
for(k in 1:K){
    fit          <- lm(V[,k] ~ g + as.matrix(pcs)[,1:10])
    mvdat[1:2,k] <- summary(fit)$coefficients[1,1:2]
    mvdat[3:4,k] <- summary(fit)$coefficients[2,1:2]
    pval[k]      <- summary(fit)$coefficients[2,4]
}
mvdat[2,]  <- mvdat[2,]^2##square the se of the intercept
mvdat[4,]  <- mvdat[4,]^2###square the se of the coefficient

print("compute effect")

DUt        <- D[1:K] * t(U[,1:K])
res.effect <- get.effectonWCs.fromPCA(mvdat=mvdat, DUt=DUt, Xmean = alpha.Mean, tcut=tcut)

mean.alpha <- res.effect$mean.alpha
var.alpha  <- res.effect$var.alpha
mean.beta  <- res.effect$mean.beta
var.beta   <- res.effect$var.beta

alpha.m <- matrix(mean.alpha, nrow=J, byrow=FALSE)
alpha.v <- matrix(var.alpha, nrow=J, byrow=FALSE)
beta.m  <- matrix(mean.beta, nrow=J, byrow=FALSE)
beta.v  <- matrix(var.beta, nrow=J, byrow=FALSE)

fit           <- lm(mu.mean~g+pcs)
total.alpha.m <- summary(fit)$coefficients[1,1]
total.alpha.v <- summary(fit)$coefficients[1,2]^2
total.beta.m  <- summary(fit)$coefficients[2,1]
total.beta.v  <- summary(fit)$coefficients[2,2]^2
pval[K+1]     <- summary(fit)$coefficients[2,4]

#res.effect.data<-reverseWTeffect.PCA(alpha.m,alpha.v,beta.m,beta.v,total.alpha.m,total.alpha.v,total.beta.m,total.beta.v)

res.effect.data <- reverseWTeffect.PCA(alpha.m,alpha.v,beta.m,beta.v,0,0,0,0)
effect.mean     <- t(as.matrix(res.effect.data$effect.mean))
effect.var      <- t(as.matrix(res.effect.data$effect.var))


pdf(file=file.path(dir.name,"effectmean.zerototalintensity.nocs.pcs.pdf"))

fra         <- 2
effect.mean <- reconstruct.1D.Matrix(effect.mean, N.0, subset.col, bin.size)
effect.var  <- reconstruct.1D.Matrix(effect.var, N.0, subset.col, bin.size)
ybottom     <- effect.mean - fra*sqrt(effect.var)
ytop        <- effect.mean + fra*sqrt(effect.var)
ymax        <- max(ytop) + 0.0000000001
ymin        <- min(ybottom) - 0.0000000001
wh.bottom   <- which(ybottom > 0)
wh.top      <- which(ytop < 0)
N.1         <- length(effect.mean)
high.wh     <- sort(unique(union(wh.bottom, wh.top)))
xval        <- 1:N.1
xmin        <- min(xval)
xmax        <- max(xval)
col.posi    <- xval[high.wh]


par(mfrow=c(3,1))
plot(ybottom, pch=".", ylim=c(ymin,ymax), main="EffectMeanofPCsonWC", col="green")
points(ytop, pch=".", col="green")
points(effect.mean, pch=".", col="dark green")
abline(h=0, col="red")

N.polygons  <- length(col.posi)
if (N.polygons > 0)
    for(j in 1:N.polygons)
        polygon(c(polygon(c(col.posi[j]-0.5, col.posi[j]-0.5, col.posi[j]+0.5, col.posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col=rgb(1, 0, 0,0.5), border = NA)))
ymin <- min(eigendeviations)
ymax <- max(eigendeviations)
plot(reconstruct.1D.Matrix(t(as.matrix(eigendeviations[,1])), N.0, subset.col, bin.size), ylim=c(ymin, ymax), ylab="y", xlab=NULL, pch=".", main="eigendeviation 1 in original space")
plotGenePred(Transcripts, locus.start, locus.end)
if (data=="YUR")
    plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)
dev.off()

##########################################
print("compute mvdat,tcut=0,totaleffectmatters")
tcut       <- 0
pval   <- rep(NA, K+1)
mvdat  <- matrix(NA, nr = 4, nc = K)
for(k in 1:K){
    fit          <- lm(V[,k] ~ g + as.matrix(pcs)[,1:10])
    mvdat[1:2,k] <- summary(fit)$coefficients[1,1:2]
    mvdat[3:4,k] <- summary(fit)$coefficients[2,1:2]
    pval[k]      <- summary(fit)$coefficients[2,4]
}
mvdat[2,]  <- mvdat[2,]^2##square the se of the intercept
mvdat[4,]  <- mvdat[4,]^2###square the se of the coefficient

print("compute effect")

DUt        <- D[1:K] * t(U[,1:K])
res.effect <- get.effectonWCs.fromPCA(mvdat=mvdat, DUt=DUt, Xmean = alpha.Mean, tcut=tcut)

mean.alpha <- res.effect$mean.alpha
var.alpha  <- res.effect$var.alpha
mean.beta  <- res.effect$mean.beta
var.beta   <- res.effect$var.beta

alpha.m <- matrix(mean.alpha, nrow=J, byrow=FALSE)
alpha.v <- matrix(var.alpha, nrow=J, byrow=FALSE)
beta.m  <- matrix(mean.beta, nrow=J, byrow=FALSE)
beta.v  <- matrix(var.beta, nrow=J, byrow=FALSE)

fit           <- lm(mu.mean~g+pcs)
total.alpha.m <- summary(fit)$coefficients[1,1]
total.alpha.v <- summary(fit)$coefficients[1,2]^2
total.beta.m  <- summary(fit)$coefficients[2,1]
total.beta.v  <- summary(fit)$coefficients[2,2]^2
pval[K+1]     <- summary(fit)$coefficients[2,4]

res.effect.data<-reverseWTeffect.PCA(alpha.m,alpha.v,beta.m,beta.v,total.alpha.m,total.alpha.v,total.beta.m,total.beta.v)

##res.effect.data <- reverseWTeffect.PCA(alpha.m,alpha.v,beta.m,beta.v,0,0,0,0)
effect.mean     <- t(as.matrix(res.effect.data$effect.mean))
effect.var      <- t(as.matrix(res.effect.data$effect.var))

pdf(file=file.path(dir.name, "effectmean.totalintensity.nocs.pcs.pdf"))
fra         <- 2
effect.mean <- reconstruct.1D.Matrix(effect.mean, N.0, subset.col, bin.size)

effect.var  <- reconstruct.1D.Matrix(effect.var, N.0, subset.col, bin.size)

ybottom     <- effect.mean - fra*sqrt(effect.var)
ytop        <- effect.mean + fra*sqrt(effect.var)
ymax        <- max(ytop) + 0.0000000001
ymin        <- min(ybottom) - 0.0000000001
wh.bottom   <- which(ybottom > 0)
wh.top      <- which(ytop < 0)
N.1         <- length(effect.mean)
high.wh     <- sort(unique(union(wh.bottom, wh.top)))
xval        <- 1:N.1
xmin        <- min(xval)
xmax        <- max(xval)
col.posi    <- xval[high.wh]


par(mfrow=c(3,1))
plot(ybottom, pch=".", ylim=c(ymin,ymax), main="EffectMeanofPCsonWC", col="green")
points(ytop, pch=".", col="green")
points(effect.mean, pch=".", col="dark green")
abline(h=0, col="red")

N.polygons  <- length(col.posi)
if (N.polygons > 0)
    for(j in 1:N.polygons)
        polygon(c(polygon(c(col.posi[j]-0.5, col.posi[j]-0.5, col.posi[j]+0.5, col.posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col=rgb(1, 0, 0,0.5), border = NA)))
ymin <- min(eigendeviations)
ymax <- max(eigendeviations)
plot(reconstruct.1D.Matrix(t(as.matrix(eigendeviations[,1])), N.0, subset.col, bin.size), ylim=c(ymin, ymax), ylab="y", xlab=NULL, pch=".", main="eigendeviation 1 in original space")
plotGenePred(Transcripts, locus.start, locus.end)
if (data=="YUR")
    plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)
dev.off()


save(list=ls(),file=file.path(dir.name, paste0(locus.name,".Robj")))

##############################################################


##################################################################

zeromat=alpha.mean[,g==0]
##zeromat=exp(zeromat)/(1+exp(zeromat))
test=apply(zeromat,1,function(x){mean(x)})
estmat=reverseWT.PCA(matrix(test, nrow=J, byrow=FALSE), threshold=NULL, highlevelfree=2, logscale=FALSE)$est
estmat <- reconstruct.1D.Matrix(t(as.matrix(estmat)), N.0, subset.col, bin.size)

onemat=alpha.mean[,g==1]
#onemat=exp(onemat)/(1+exp(onemat))
test1=apply(onemat,1,function(x){mean(x)})
estmat1=reverseWT.PCA(matrix(test1, nrow=J, byrow=FALSE), threshold=NULL, highlevelfree=2, logscale=FALSE)$est
estmat1 <- reconstruct.1D.Matrix(t(as.matrix(estmat1)), N.0, subset.col, bin.size)



twomat=alpha.mean[,g==2]
#twomat=exp(twomat)/(1+exp(twomat))	
test2=apply(twomat,1,function(x){mean(x)})
estmat2=reverseWT.PCA(matrix(test2, nrow=J, byrow=FALSE), threshold=NULL, highlevelfree=2, logscale=FALSE)$est
estmat2 <- reconstruct.1D.Matrix(t(as.matrix(estmat2)), N.0, subset.col, bin.size)

pdf(file=file.path(dir.name,"retransformedbygenotype.pdf"))
par(mfrow=c(4,1))
plot(estmat,col="red",cex=0.5)#,ylim=c(0,4e-4))
plot(estmat1,col="blue",cex=0.5)#,ylim=c(0,4e-4))
plot(estmat2,col="green",cex=0.5)#,ylim=c(0,4e-4))
plotGenePred(Transcripts, locus.start, locus.end)
if (data=="YUR")
    plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)
plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)
dev.off()

zeromat=M.0[g==0,]
test=apply(zeromat,2,function(x){mean(x)})
onemat=M.0[g==1,]
test1=apply(onemat,2,function(x){mean(x)})
twomat=M.0[g==2,]
test2=apply(twomat,2,function(x){mean(x)})

ymin=min(c(test,test1,test2))
ymax=max(c(test,test1,test2))

pdf(file=file.path(dir.name,"originaldata.pdf"))
par(mfrow=c(3,1))
plot(test,col="red",cex=0.5,ylim=c(ymin,ymax))
plot(test1,col="pink",cex=0.5,ylim=c(ymin,ymax))
plot(test2,col="yellow",cex=0.5,ylim=c(ymin,ymax))
dev.off()

save(list=ls(),file=file.path(dir.name, paste0(locus.name,".Robj")))
############


res.est.PB = cyclespin.smooth(M$M,g=g, fast.approx=TRUE, read.depth = read.depth,lm.approx=TRUE)


effect.cs.mean <- t(as.matrix(res.est.PB$effect.mean))
effect.cs.var <- t(as.matrix(res.est.PB$effect.var))

fra         <- 2
effect.cs.mean <- reconstruct.1D.Matrix(effect.cs.mean, N.0, subset.col, bin.size)
effect.cs.var <- reconstruct.1D.Matrix(effect.cs.var, N.0, subset.col, bin.size)
fra = 2
beta_l = effect.cs.mean - fra*sqrt(effect.cs.var)
beta_r = effect.cs.mean + fra*sqrt(effect.cs.var)

ymin_beta = min(beta_l) - 0.0000000001
ymax_beta = max(beta_r) + 0.0000000001

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

col_posi = xval[high_wh]


pdf(file=file.path(dir.name,"cyclespin.pdf"))
par(mfrow=c(3,1))
plot(beta_l, pch=".", ylim=c(ymin_beta,ymax_beta), main="EffectMeanofCSonWC", col="green")
points(beta_r,pch=".", col="green")
points(effect.cs.mean, pch=".", col="dark green")
abline(h=0, col="red")

N.polygons  <- length(col_posi)
if (N.polygons > 0)
    for(j in 1:N.polygons)
        polygon(c(polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col=rgb(1, 0, 0,0.5), border = NA)))
ymin <- min(eigendeviations)
ymax <- max(eigendeviations)
plot(reconstruct.1D.Matrix(t(as.matrix(eigendeviations[,1])), N.0, subset.col, bin.size), ylim=c(ymin, ymax), ylab="y", xlab=NULL, pch=".", main="eigendeviation 1 in original space")
plotGenePred(Transcripts, locus.start, locus.end)
if (data=="YUR")
    plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)
dev.off()








###Try to confirm with cyclespin estimate of effect size

#read.depth=read.depth[which(!is.na(g))]


##res = cyclespin.smooth(M2,g=g)
#res2 = cyclespin.smooth(M$M,g, read.depth = read.depth)

#pdf(file=file.path(dir.name, "effectmean.wCS.nopolygon.pdf"))
#par(mfrow=c(4,1))

#par(mar = c(2,4,2,0))
#plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin, ymax),xlim =c(xmin, xmax), main ="Without cyclespin", axes=FALSE)
#axis(2)
#N.polygons  <- length(col.posi)
#if (N.polygons > 0)
#    for(j in 1:N.polygons)
#        polygon(c(polygon(c(col.posi[j]-0.5, col.posi[j]-0.5, col.posi[j]+0.5, col.posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col=rgb(1, 0, 0,0.5), border = #NA)))
#points (1:N.1, ybottom, type="l", ylim=c(ymin,ymax), main="EffectMeanofPCsonWC", col="green")
#abline(h=0,col="red")

#points(1:N.1, ytop, type="l", col="green")
#points(1:N.1, effect.mean, type="l", col="dark green")
#plot(1:N.1,M.PC1,type="l")
#plotGenePred(Transcripts, locus.start, locus.end)
#if (data=="YUR")
 #   plotListLociLine(locus.start, locus.end, list.loci.file, list.loci.line)


#effect.mean <- reconstruct.1D.Matrix(t(as.matrix(res2$effect.mean)), N.0, subset.col, filtering)
#effect.var  <- reconstruct.1D.Matrix(t(as.matrix(res2$effect.var)), N.0, subset.col, filtering)
#ybottom     <- effect.mean - fra*sqrt(effect.var)
#ytop        <- effect.mean + fra*sqrt(effect.var)
#ymax        <- max(ytop) + 0.0000000001
#ymin        <- min(ybottom) - 0.0000000001
#wh.bottom   <- which(ybottom > 0)
#wh.top      <- which(ytop < 0)
#N.1         <- length(effect.mean)


#high.wh     <- sort(unique(union(wh.bottom, wh.top)))
#xval        <- 1:N.1
#xmin        <- min(xval)
#xmax        <- max(xval)
#col.posi    <- xval[high.wh]

#plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin, ymax),xlim =c(xmin, xmax), main ="With_cyclespin", axes=FALSE)

#N.polygons  <- length(col.posi)
#if (N.polygons > 0)
#    for(j in 1:N.polygons)
#        polygon(c(polygon(c(col.posi[j]-0.5, col.posi[j]-0.5, col.posi[j]+0.5, col.posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col=rgb(1, 0, 0,0.5), border = NA)))
#abline(h=0,col="red")
#points(1:N.1, ybottom, type="l", ylim=c(ymin,ymax), main="EffectMeanCConWC", col="green")
#points(1:N.1, ytop, type="l", col="green")
#points(1:N.1, effect.mean, type="l", col="dark green")

#dev.off()

save(list=ls(),file=file.path(dir.name, paste0(locus.name,".Robj")))


