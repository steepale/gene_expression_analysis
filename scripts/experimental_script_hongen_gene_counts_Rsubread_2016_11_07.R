
# Bam file input: 017834-2_paired_norRNA_sorted.bam

biocLite("Rsubread")
library("Rsubread")



fc <- featureCounts("834.bam", annot.ext="Gallus_gallus.Gallus_gallus-5.0.86.gtf", isGTFAnnotationFile=TRUE, GTF.attrType="gene_id")
names(fc)
colnames(fc$counts) <- rownames(targets)

classifications <- read.delim("targets_failures.txt", stringsAsFactors=FALSE)

classifications

samples <- c("017834-2")

GCs <- read.delim("017834-2_gene_counts2.txt", row.names="ENSEMBL")


read.sample <- function(sample.name) {
        file.name <- paste(sample.name, "_gene_counts2.txt", sep="")
        result <- read.delim(file.name, col.names=c("SYMBOL", "count"), sep="\t", colClasses=c("character", "numeric"), row.names=1)
}

sample.1 <- read.sample(samples[1])

classifications <- read.delim("targets_failures.txt", stringsAsFactors=FALSE)

classifications

group <- paste(classifications$group, classifications$description, sep=".")
group <- factor(group)
table(group)








####
obj_test <- getBM(attributes=c('entrezgene','external_gene_name'), 
                  filters = 'external_gene_name', 
                  values = rownames(all.data), mart = ensembl)

for (i in 2:length(rownames(tr))){
        print(i)
        gene<-rownames(tr)[i]
        rownames(tr)[i]<-obj_test[obj_test$external_gene_name==gene,]$external_gene_name
}

####


