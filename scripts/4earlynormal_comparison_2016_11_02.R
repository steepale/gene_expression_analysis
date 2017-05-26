targets <- read.delim("targets_4earlynormal.txt", stringsAsFactors=FALSE)

targets

group <- paste(targets$group, targets$description, sep=".")
group <- factor(group)
table(group)

samples <- c("017733", "017748", "017820", "017824")

read.sample <- function(sample.name) {
        file.name <- paste(sample.name, "_gene_counts.txt", sep="")
        result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"), row.names=1)
}

sample.1 <- read.sample(samples[1])
sample.2 <- read.sample(samples[2])
sample.3 <- read.sample(samples[3])
sample.4 <- read.sample(samples[4])

all.data <- data.frame(sample.1, sample.2$count, sample.3$count, sample.4$count)

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

AvsB <- makeContrasts(NCA.NormalA-NCB.NormalB, levels=design)

res <- glmQLFTest(fit, contrast=AvsB)

is.de <- decideTestsDGE(res)

summary(is.de)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")


tr <- glmTreat(fit, contrast=NvsT, lfc=log2(1.5))
topTags(tr, n=27)

is.de <- decideTestsDGE(tr)
summary(is.de)

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")

logCPM <- cpm(dge, prior.count=2, log=TRUE)
rownames(logCPM) <- dge$genes$Symbol
colnames(logCPM) <- paste(dge$samples$group, 1:2, sep="-")

o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:30],]

logCPM <- t(scale(t(logCPM)))

library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))


query = perform.query(c("LOC107055483", "LOC374195", "LOC100859467", "UBAP2", "LOC427025", "LOC107049271", "LOC107049285", "LOC107049047", "LOC107049297", "LOC107049279", "LOC107049174", "LOC107052259", "UBE2R2L", "LOC100859602", "LOC107049296", "LOC426615", "LOC100859779", "LOC107049306", "LOC107052275", "LOC100857891", "LOC101752211", "LOC430766", "LOC107049287", "LOC768814", "LOC431003", "SUB1L2", "LOC426344"))



