library(DESeq2)
library(edgeR)
library(sva)

#biocLite("bladderbatch")

# get expression matrix
samples <- read.delim("/mnt/projects/iamp/results/anduril/execute/_qcReport_counts_array1/array/_index")
rownames(samples) <- samples$Key

samples$type <- NA
samples$type[grepl("^C\\d", rownames(samples))] <- "ER"
samples$type[grepl("^A\\d", rownames(samples))] <- "iAMP21"
samples$type[grepl("^D\\d", rownames(samples))] <- "DS"
samples$type[grepl("^S\\d", rownames(samples))] <- "CD19+"
samples$type <- as.factor(samples$type)

cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="/", design=~1)
expressed <- rowSums(counts(cds)) >= 10
dge <- DGEList(counts=counts(cds[expressed,]))
dge.norm <- calcNormFactors(dge, method="TMM")
y <- voom(dge.norm)

# get confounding variable
lib <- read.delim("/mnt/projects/iamp/results/qc/allpatients.stats.txt")
rownames(lib) <- sub(".*C57C3ACXX_CV_([^_]+)_.*", "\\1", lib$sample)
lib$pct_rrna <- (lib$exonic-lib$non.rRNA)/lib$exonic
lib$pct_rna <- (lib$non.rRNA-lib$protein)/lib$non.rRNA
lib$pct_intronic <- (lib$mapped-lib$exonic)/lib$mapped
lib$pct_dup <- (lib$uniquely.mapped-lib$non.duplicates)/lib$uniquely.mapped
pca <- prcomp(lib[,c("pct_rrna", "pct_rna", "pct_intronic", "pct_dup")], center = TRUE, scale. = TRUE) 
lib$PC1 <- pca$x[,1]
lib$PC2 <- pca$x[,2]

des <- merge(samples[,c("Key", "type")], lib[,c("PC1", "PC2")], by.x="Key", by.y="row.names")

# for control experiment
#des$PC1 <- gtools::permute(des$PC1)
#des$PC2 <- gtools::permute(des$PC2)

des$batchPC1 <- NA
des$batchPC1[des$PC1 < mean(des$PC1)-sd(des$PC1)] <- "low"
des$batchPC1[des$PC1 > mean(des$PC1)+sd(des$PC1)] <- "high"
des$batchPC1[is.na(des$batchPC1)] <- "norm"

des$batchPC2 <- NA
des$batchPC2[des$PC2 < mean(des$PC2)-sd(des$PC2)] <- "low"
des$batchPC2[des$PC2 > mean(des$PC2)+sd(des$PC2)] <- "high"
des$batchPC2[is.na(des$batchPC2)] <- "norm"

#1 - correct for PC1 
mod1 <- model.matrix(~type+batchPC2, data=des) 
bat1 <- ComBat(y$E, des$batchPC1, mod1)

#1 - correct for PC2
mod2 <- model.matrix(~type, data=des) 
bat2 <- ComBat(bat1, des$batchPC2, mod2)


# plot heatmap
library(gplots)
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
dists <- as.matrix(as.dist(1-cor(bat2, method="spearman"))) # spearman correlation
pdf("/mnt/projects/iamp/results/heatmap-voom-batchcorrected.pdf")
heatmap.2(dists, trace="none", scale="none", col=rev(hmcol), margin=c(7, 7), cexRow=0.8, cexCol=0.8, key.title="")
dev.off()

