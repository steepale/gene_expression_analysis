

targets <- read.delim("targets_controls.txt", stringsAsFactors=FALSE)

targets

group <- paste(targets$group, targets$description, sep=".")
group <- factor(group)
table(group)

samples_all <- c("017733", "017748", "017820", "017824", "017936", "017939", "017945", "017947")

read.sample <- function(sample.name) {
        file.name <- paste(sample.name, "_gene_counts.txt", sep="")
        result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"), row.names=1)
}

sample.1 <- read.sample(samples[1])
sample.2 <- read.sample(samples[2])
sample.3 <- read.sample(samples[3])
sample.4 <- read.sample(samples[4])
sample.5 <- read.sample(samples[5])
sample.6 <- read.sample(samples[6])
sample.7 <- read.sample(samples[7])
sample.8 <- read.sample(samples[8])


all.data <- data.frame(sample.1, sample.2$count, sample.3$count, sample.4$count, sample.5$count, sample.6$count, sample.7$count, sample.8$count)


colnames(all.data)[1:ncol(all.data)] <- samples

dge = DGEList(counts=all.data, group=group)
dge$samples

keep <- rowSums(cpm(dge) > 0.5) >= 2
table(keep)

dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)
dge$samples


pch <- c(0,1)
plotMDS(dge, pch=pch[group])
legend("topleft", legend=levels(group), pch=pch)

plotMD(dge, column=1)
abline(h=0, col="red", lty=2, lwd=2)

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

dge <- estimateDisp(dge, design, robust=TRUE)

plotBCV(dge)

fit <- glmQLFit(dge, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)

EvsL <- makeContrasts(E.Early-L.Late, levels=design)

res <- glmQLFTest(fit, contrast=EvsL)

is.de <- decideTestsDGE(res)

summary(is.de)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")


tr <- glmTreat(fit, contrast=EvsL, lfc=log2(1.5))
topTags(tr, n=40)

is.de <- decideTestsDGE(tr)
summary(is.de)

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")

logCPM <- cpm(dge, prior.count=2, log=TRUE)
rownames(logCPM) <- dge$genes$Symbol
colnames(logCPM) <- paste(dge$samples$group, 1:2, sep="-")

o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:40],]

logCPM <- t(scale(t(logCPM)))

library(gplots)
col.pan <- colorpanel(40, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))



query = perform.query(c("LOC107057378", "TARP", "CTSG", "LOC423462", "XCL1", "IL2RA", "UBTD2", "LOC396098", "CSF2RB", "RP11-47I22.4", "PLK2", "LOC101748259", "LOC101750578", "IL2RB", "HSP90AA1", "DNAJA4", "SOX13", "ST18", "FREM3", "LIN28B", "C3H8ORF80", "GRIN3A", "CCL20", "FREM2", "CADPS", "GPR82", "ZNF502", "INHBA", "ATP10A", "LOC107054604", "HDC", "OTUD1", "ARHGEF9", "SLIT3", "B3GNT4", "THBS2", "ZNF366", "ABI3BP", "HSPH1", "LOC107053544"))

