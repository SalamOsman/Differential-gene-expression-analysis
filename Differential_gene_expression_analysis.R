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
library(factoextra)
library(viridis)
library(EnhancedVolcano)

# Importing the count data
countData <- read.table('C:/Users/khan_/Desktop/DGEA/HNSC.mRNAseq_raw_counts.txt', 
                        header = TRUE, sep = "\t", check.names = FALSE)
head(countData[1:10])

# Renaming and reshaping of the IDs column in countData.
colnames(countData)[1] <- "ID"
IDs <- data.frame(do.call("rbind", strsplit(as.character(countData$ID), "|", fixed = TRUE)))
IDs$X2 <- NULL
colnames(IDs)[1] <- "Gene_symbol"

# Building a countData dataframe with Gene_symbol as IDs by filtering. First, let's reshape the countData accordingly to the matching IDs.
countData <- cbind(countData, IDs)
countData <- countData[countData$Gene_symbol != "?", ]
countData <- countData[!duplicated(countData$Gene_symbol), ]
rownames(countData) <- countData$Gene_symbol

# Discarding the unwanted column from the expression matrix.
countData[ ,c('ID', 'Gene_symbol')] <- list(NULL)
head(countData[1:10])

# Importing the metafile with clinical information
metaData <- read.table('C:/Users/khan_/Desktop/DGEA/meta_data.txt', 
                       header = TRUE, sep = "\t")
head(metaData)

# Resetting metaData file containing clinical information for the gene expression matrix (countData). 
countData <- as.data.frame(t(countData))
countData <- filter(countData , row.names(countData) %in% metaData$IDs)
countData <- as.data.frame(t(countData))

# Checking the library sizes of the samples in the raw counts file.
librarySizes <- colSums(countData)
barplot(librarySizes, las=2, border = "NA", main="Library sizes of cancer vs normal samples")

# Checking the distribution of the counts in samples. 
logcounts <- log2(countData + 1)
statusCol <- as.numeric(factor(metaData$Sample_type) + 1)
boxplot(logcounts, xlab="", ylab="Log2(Counts)", las=2,col=statusCol, main = "raw Count distribution in cancer vs normal samples")
abline(h=median(as.matrix(logcounts)), col="red")

# PCA of the samples based on sample type.
pcDat <- prcomp(t(countData), center = T)
fviz_pca_ind(pcDat, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = metaData$Sample_type, col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black", repel = TRUE, legend.title = "Sample types") + 
  ggtitle("Cancer vs Normal samples based on types") + 
  theme(plot.title = element_text(hjust = 0.5))

# PCA of the samples based on cancer stages
fviz_pca_ind(pcDat, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = metaData$Stage, col.ind = "black", palette = "jco", addEllipses = TRUE, label = "var", col.var = "black", repel = TRUE, legend.title = "Sample types") + 
  ggtitle("Cancer vs Normal samples based on cancer stages") + 
  theme(plot.title = element_text(hjust = 0.5))

# Verifying the sample names from meta data in raw counts file.
all(metaData$IDs == colnames(countData))

# Building study designs for downstream analysis.
design <- as.formula(~ Sample_type)

# Creating dds object for each study separately and normalization.
dds <- DESeqDataSetFromMatrix(countData =  round(countData), colData = metaData, design = design)

# Dropping/discarding low counts genes (< 5) from the studies.
keep <- rowSums(counts(dds)) > 5
dds <- dds[keep,]
nrow(dds)

# Calculating distance among the samples.
vsd <- vst(dds, blind = FALSE)
colData(vsd)
sampleDists <- dist(t(assay(vsd)))

# Plotting a heatmap of the distance score
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$Status, vsd$Age, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists, show_rownames = F, show_colnames = F,
         clustering_distance_cols = sampleDists, cluster_col = F, cluster_rows = F,
         col = magma(10), main = "Cancer vs Control")

# Post-normalization of the dds object.
dds <- estimateSizeFactors(dds)

# Performing differential gene expression analysis using DESEQ function.
dds <- DESeq(dds)
res <- results(dds)
res

# Making a contrast group and summarizing the results.
res_group <- results(dds, contrast=c("Sample_type", "HNSCC", "Normal"))
mcols(res_group, use.names = TRUE)
summary(res_group)

# Plotting volcano charts of the studies. 
EnhancedVolcano(res_group, lab = rownames(res), x = 'log2FoldChange', y = 'padj', 
                FCcutoff = 2, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','log2FC','adj.P',
                               'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, xlim = c(-10,10), ylim = c(0, 15),
                title = "Cancer vs Normal samples", subtitle = NULL)

# Summarizing the results and sorting on the basis of expression values and significance (adjusted p value <0.05 and absolute log 2 fold change greater than 2). 
resSig <- subset(res_group, padj < 0.05 & abs(log2FoldChange) > 2)
dim(resSig)

# Shrinking the log2FC values and constructing MA plot.
resultsNames(dds)
res <- lfcShrink(dds, coef="Sample_type_Normal_vs_HNSCC", type="apeglm")
dim(res)
plotMA(res, ylim = c(-5, 5))

# Writing the list of DEGs. 
write.csv(resSig, file = "C:/Users/khan_/Desktop/DGEA/DGE_list.csv", quote = F)
