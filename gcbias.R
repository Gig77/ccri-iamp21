options(warn=1)
library(RColorBrewer)
library(DESeq2)

# get read counts
files <- list.files(path="/mnt/projects/iamp/results/htseq/", pattern=".count$")
names <- sub(".*CV_([^_]+).*", "\\1", files)
samples <- data.frame(name=names, file=files, stringsAsFactors=F)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory="/mnt/projects/iamp/results/htseq/", design=~1)
cds <- estimateSizeFactors(cds)
counts <- round(counts(cds, normalized=T))

# batches
lane <- as.factor(paste("lane", sub(".*_lane(\\d).*", "\\1", files)))
plex <- ifelse(names %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7"), "plex 1", NA)
plex <- ifelse(names %in% c("C8", "C9", "C10", "C11", "C12", "A1", "A2"), "plex 2", plex)
plex <- ifelse(names %in% c("A3", "A4", "A5", "A7", "A8", "A9", "A10"), "plex 3", plex)
plex <- ifelse(names %in% c("A11", "A12", "A13", "D1", "D2", "D3", "D4"), "plex 4", plex)
plex <- ifelse(names %in% c("D5", "D6", "D7", "D8", "S1", "S2", "S3"), "plex 5", plex)
plex <- as.factor(plex)

# filter for expressed transcripts
counts.expressed <- counts[rowMeans(counts)>10,]
#counts.expressed <- counts[apply(counts, 1, min)>10,]
nrow(counts.expressed)

# get GC content and length of (longest) transcript per gene
library("biomaRt")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75
result <- select(mart, keys=rownames(counts), keytype="ensembl_gene_id", columns=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"))
result.longest <- result[ave(result$transcript_length, result$ensembl_gene_id, FUN=max) == result$transcript_length,]
feature <- data.frame(row.names=result.longest$ensembl_gene_id, gc=result.longest$percentage_gc_content, length=result.longest$transcript_length)
nrow(feature)

# get genes that are expressed AND for which we have GC content/transcript length information
common <- intersect(rownames(counts.expressed), rownames(feature))
length(common)

library(EDASeq)
data <- newSeqExpressionSet(counts=counts.expressed[common,], featureData=feature[common,], phenoData=data.frame(lane=lane, plex=plex, row.names=colnames(counts.expressed)))

pdf("/mnt/projects/iamp/results/gcbias.pdf", width=10, height=15)
par(mfrow=c(2,1))

# plot normalized expression by sample
colors.lane <- brewer.pal(length(levels(lane)), "Accent")[lane]
colors.plex <- c("darkgray", "red", "blue", "orange", "lightgray")[plex]
boxplot(data, col=colors.lane, las=2)
boxplot(data, col=colors.plex, las=2)

# gc plot lane 1 vs. lane 2 (both mostly iAMP); variable but no systematic difference by lane
biasPlot(data[,lane %in% c("lane 1", "lane 2")], "gc", log=TRUE, ylim=c(3.5, 5.5), color_code=1, main="iAMP21")

# gc plot lane 3 vs. lane 4 (both mostly ER); variable but no systematic difference by lane
biasPlot(data[,lane %in% c("lane 3", "lane 4")], "gc", log=TRUE, ylim=c(3.5, 5.5), color_code=1, main="ER")

# gc plot lane 1 vs. lane 4 (iAMP vs. ER); systematic difference, but could be biological not technical
biasPlot(data[,lane %in% c("lane 1", "lane 4")], "gc", log=TRUE, ylim=c(3.5, 5.5), color_code=1, main="iAMP vs. ER")

# gc plot plex 1 vs. plex 2 (both mostly ER); highly variable, but mixed
biasPlot(data[,plex %in% c("plex 1", "plex 2")], "gc", log=TRUE, ylim=c(3.5, 5.5), color_code=2, main="ER")

# gc plot plex 3 vs. plex 1 (iAMP vs. ER); systematic difference, but could be biological not technical
biasPlot(data[,plex %in% c("plex 3", "plex 1")], "gc", log=TRUE, ylim=c(3.5, 5.5), color_code=2, main="iAMP vs. ER")

# bias plots for two selected samples with very different GC bias
biasPlot(data[,c("A11", "C5")], "gc", log=TRUE, ylim=c(3.5, 5.5), color_code=1, main="A11 vs. C5")
biasBoxplot(log(counts(data)[,"A11"]+0.1) - log(counts(data)[,"C5"]+0.1), fData(data)$gc, ylim=c(-5, 5), main="A11 vs. C5")

dev.off()
