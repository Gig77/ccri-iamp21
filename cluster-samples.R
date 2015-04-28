options(warn=1)
library("DESeq2")
library("RColorBrewer")
library("ggplot2")
library("gplots")
library("edgeR")

plot_heatmap <- function(data, file) {
	distsRL <- dist(t(assay(data)))
	mat <- as.matrix(distsRL)
	hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
	pdf(file, width=12, height=12)
	heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
	dev.off()
} 

plot_pca <- function(data, file) {
	pdf(file)
	p <- plotPCA(data, intgroup=c("type"), returnData=TRUE)
	percentVar <- round(100 * attr(p, "percentVar"))
	p <- ggplot(p, aes(PC1, PC2, color=type)) +
			theme_bw() +
			geom_point(size=5) + 
			geom_text(aes(label=names), size=1.5, colour="white", show_guide=F) +
			scale_colour_manual(values=c("red", "blue", "black", "orange")) +
			xlab(paste0("PC1: ", percentVar[1], "% variance")) +
			ylab(paste0("PC2: ", percentVar[2], "% variance"))
	print(p)	
	dev.off()
} 

# get counts
#------------
files <- list.files(path="~/iamp/results/htseq/", pattern=".count$")
samples <- data.frame(name=gsub(".*CV_(.\\d+)_.*", "\\1", files), file=files, type=as.character(NA), stringsAsFactors=F)
samples$type[grepl("A\\d", samples$name)] <- "iAMP" 
samples$type[grepl("C\\d", samples$name)] <- "E/R" 
samples$type[grepl("D\\d", samples$name)] <- "P/C" 
samples$type[grepl("S\\d", samples$name)] <- "CD19+" 
samples$type <- as.factor(samples$type)

cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="~/iamp/results/htseq", design=~1)
#cds <- cds[rowSums(cpm(DGEList(counts=counts(cds))) > 5 ) >= 2,] # keep only expressed genes 

#---------------------
# heatmap and pca of all expressed genes
#---------------------
rld <- rlog(cds)
plot_heatmap(rld, "~/iamp/results/sample-dist.heatmap.deseq2.rlog.pdf")
plot_pca(rld, "~/iamp/results/sample-dist.pca.deseq2.rlog.pdf")

#pdf("~/fikret/results/sample-dist.heatmap.deseq2.rlog.customclustering.pdf", width=12, height=12)
#heatmap.2(as.matrix(distsRL), hclustfun=function(x) hclust(x, method="average"), distfun=function(x) as.dist(x), trace="none", col=rev(hmcol), margin=c(13, 13))
#dev.off()

#---------------------
# heatmap and pca of expressed chr21 genes
#---------------------
load("~/generic/data/ensembl/genes.GRCh37v75.biomart.RData")
genes <- genes[genes$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"),]
transcripts <- data.frame(ensembl_gene_id=rownames(counts(cds)))
transcripts <- merge(transcripts, genes[,c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "start_position", "end_position")], by="ensembl_gene_id", all.x=T) # add gene annotation
transcripts <- transcripts[!duplicated(transcripts$ensembl_gene_id),]

rld.21 <- rlog(cds[!is.na(transcripts$chromosome_name) & transcripts$chromosome_name=="21",])
plot_heatmap(rld.21, "~/iamp/results/sample-dist.heatmap.chr21.deseq2.rlog.pdf")
plot_pca(rld.21, "~/iamp/results/sample-dist.pca.chr21.deseq2.rlog.pdf")

#---------------------
# heatmap for miR target genes
#---------------------
load("~/generic/data/ensembl/genes.GRCh37v75.biomart.RData")
genes <- genes[genes$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"),]
ens2sym <- data.frame(ensembl_gene_id=rownames(cds))
ens2sym <- merge(ens2sym, genes[,c("ensembl_gene_id", "hgnc_symbol")], all.x=T)
ens2sym <- ens2sym[!duplicated(ens2sym$ensembl_gene_id),]

# validated targets
validated <- read.delim("~/iamp/results/miRecords.v4.validatedTargets.txt")
validated.hgnc <- as.character(unique(validated[grepl("let-7c|miR-99a|miR-125b", validated$miRNA_mature_ID),"symbol"]))
validated.ensembl <- ens2sym[ens2sym$hgnc_symbol %in% validated.hgnc, "ensembl_gene_id"]

hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
mat <- assay(rld)[validated.ensembl,]
rownames(mat) <- ens2sym[match(rownames(mat), ens2sym$ensembl_gene_id), "hgnc_symbol"]
pdf("~/iamp/results/sample-dist.heatmap.let-7c.miR-99a.miR-125b.validated-targets.deseq2.rlog.pdf", width=8, height=12)
heatmap.2(mat, col=hmcol, dendrogram="both", trace="none", margin=c(10, 6))
dev.off()

# validated + predicted targets
predicted <- read.delim("~/iamp/results/miRecords.v4.predictedTargets.by5programs.txt")
predicted.hgnc <- as.character(unique(predicted[grepl("let-7c|miR-99a|miR-125b", predicted$miRNA.ID),"Symbol"]))
predicted.ensembl <- ens2sym[ens2sym$hgnc_symbol %in% predicted.hgnc, "ensembl_gene_id"]

hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
mat <- assay(rld)[unique(c(validated.ensembl, predicted.ensembl)),]
rownames(mat) <- ens2sym[match(rownames(mat), ens2sym$ensembl_gene_id), "hgnc_symbol"]
pdf("~/iamp/results/sample-dist.heatmap.let-7c.miR-99a.miR-125b.validated_plus_predicted-targets.deseq2.rlog.pdf", width=8, height=24)
heatmap.2(mat, col=hmcol, dendrogram="both", trace="none", margin=c(10, 6), cexRow=0.2)
dev.off()

#---------------------
# heatmap and pca of chr1 genes
#---------------------
load("~/generic/data/ensembl/genes.GRCh37v75.biomart.RData")
genes <- genes[genes$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"),]
transcripts <- data.frame(ensembl_gene_id=rownames(counts(cds)))
transcripts <- merge(transcripts, genes[,c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "start_position", "end_position")], by="ensembl_gene_id", all.x=T) # add gene annotation
transcripts <- transcripts[!duplicated(transcripts$ensembl_gene_id),]

rld.chr1 <- rlog(cds[!is.na(transcripts$chromosome_name) & transcripts$chromosome_name=="1",])
plot_heatmap(rld.chr1, "~/iamp/results/sample-dist.heatmap.chr1.deseq2.rlog.pdf")
plot_pca(rld.chr1, "~/iamp/results/sample-dist.pca.chr1.deseq2.rlog.pdf")
