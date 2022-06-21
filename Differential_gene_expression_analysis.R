# Loading required libraries for the analysis
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

#importing the count data
countData <- read.table('Path_to_repository/HNSC.mRNAseq_raw_counts.txt', 
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(countData)

#importing the metafile with clinical information
metaData <- read.table('Path_to_repository/HNSCC_DEGs/meta_data.txt', 
                       header = TRUE, sep = "\t", row.names = 1)
metaData

# Matching up the samples in count data with metafile
all(rownames(metaData) %in% colnames(countData))
all(rownames(metaData) == colnames(countData))
countData <- countData[, rownames(metaData)]
all(rownames(metaData) == colnames(countData))
dim(countData)
dim(metaData)

# Grouping the count data based clinical information
Status <- factor(metaData$Sample_type)
Age <- factor(metaData$Alcohol_consumption)
Gender <- factor(metaData$Gender)
Clinical <- factor(metaData$Clinical_stage)

# Create the design formula for the groups in current study
design <- as.formula(~Status)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = design)

# Filtering (We will keep all genes where the total number of reads across all samples is greater than 5.)
nrow(dds)
keep <- rowSums(counts(dds)) > 5
dds <- dds[keep,]
nrow(dds)

# Normalization of the variance in count data across samples.
vsd <- vst(dds, blind = FALSE)
colData(vsd)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Calculating pairwise distance score across the samples
sampleDists <- dist(t(assay(vsd)))
sampleDists

# Plotting a heatmap of the distance score
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

# Performing differential gene expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

# Alternative to the above chunk
res <- results(dds, contrast=c("Status","Infected","Negative"))
mcols(res, use.names = TRUE)
summary(res)

# Summarizing the results and sorting on the basis of expression values adn significance. 
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC2 <- results(dds, lfcThreshold=1)
table(resLFC2$padj < 0.1)
sum(res$padj < 0.05, na.rm=TRUE)
resSig <- subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
dim(resSig)

# Shrinking the log2FC values.
resultsNames(dds)
res <- lfcShrink(dds, coef="Status_Negative_vs_Infected", type="apeglm")
dim(res)
plotMA(res, ylim = c(-5, 5))

# Plotting the distogram of p-values
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

# Sorting the list based on P-values
resOrdered <- res[order(res$pvalue),]
dim(resOrdered)
head(resOrdered)

# Writing the list of DEGs. 
write.csv(resLFC2, file = "Path_to_repository/DGE_list.csv", quotes = F)
