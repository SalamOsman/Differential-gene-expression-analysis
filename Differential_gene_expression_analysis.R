library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(Glimma)
library(ggfortify)
library(genefilter)

#inputing data
countData <- read.table('C:/Users/Salam/Desktop/GitHub/HNSCC_DEGs/HNSC.mRNAseq_raw_counts.txt', 
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(countData)

#inputing metadata
metaData <- read.table('C:/Users/Salam/Desktop/GitHub/HNSCC_DEGs/meta_data.txt', 
                       header = TRUE, sep = "\t", row.names = 1)
metaData

#matching up
all(rownames(metaData) %in% colnames(countData))
all(rownames(metaData) == colnames(countData))
countData <- countData[, rownames(metaData)]
all(rownames(metaData) == colnames(countData))
dim(countData)
dim(metaData)

#grouping data
Status <- factor(metaData$Sample_type)
Age <- factor(metaData$Alcohol_consumption)
Gender <- factor(metaData$Gender)
Clinical <- factor(metaData$Clinical_stage)

# create the design formula for group
design <- as.formula(~Status)

# create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = design)

#filtering (We will keep all genes where the total number of reads across all samples is greater than 10.)
nrow(dds)
keep <- rowSums(counts(dds)) > 5
dds <- dds[keep,]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

#important
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#important
normalized_counts <- counts(dds, normalized=TRUE)

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$Status, vsd$Age, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

rownames(sampleDistMatrix) <- paste( vsd$Status, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

dds <- DESeq(dds)
res <- results(dds)
res

res <- results(dds, contrast=c("Status","Infected","Negative"))
mcols(res, use.names = TRUE)
summary(res)

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC2 <- results(dds, lfcThreshold=1)
table(resLFC2$padj < 0.1)

sum(res$padj < 0.05, na.rm=TRUE)

resSig <- subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
dim(resSig)

resultsNames(dds)
res <- lfcShrink(dds, coef="Status_Negative_vs_Infected", type="apeglm")
dim(res)
plotMA(res, ylim = c(-5, 5))

hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

resOrdered <- res[order(res$pvalue),]
dim(resOrdered)
head(resOrdered)




write.csv(resLFC2, file = "C:/Users/Salam/Desktop/sars/GSE152075/DGE.csv")
