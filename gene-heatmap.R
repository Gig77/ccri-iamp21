library(RColorBrewer)
library("gplots")

iAMP.vs.ER <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-ER.tsv")
iAMP.vs.DS <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-DS.tsv")
ER.vs.DS <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-DS.tsv")
iAMP.vs.immature <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-immature.tsv")
DS.vs.immature <- read.delim("/mnt/projects/iamp/results/deseq/DS-vs-immature.tsv")
ER.vs.immature <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-immature.tsv")
iAMP.vs.preB <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-preB.tsv")
DS.vs.preB <- read.delim("/mnt/projects/iamp/results/deseq/DS-vs-preB.tsv")
ER.vs.preB <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-preB.tsv")
iAMP.vs.mature <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-mature.tsv")
DS.vs.mature <- read.delim("/mnt/projects/iamp/results/deseq/DS-vs-mature.tsv")
ER.vs.mature <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-mature.tsv")

cols.keep <- c("id", "hgnc_symbol", "description", "chromosome_name", "start_position", "end_position", "log2FoldChange", "padj")
cols.merge <- c("id", "hgnc_symbol", "description", "chromosome_name", "start_position", "end_position")

merged <- merge(iAMP.vs.ER[,cols.keep], iAMP.vs.DS[,cols.keep], by=cols.merge, suffixes=c(".iAMP.vs.ER", ".iAMP.vs.DS"), all=T)
merged <- merge(merged, ER.vs.DS[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".ER.vs.DS")
merged <- merge(merged, iAMP.vs.immature[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".iAMP.vs.immature")
merged <- merge(merged, DS.vs.immature[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".DS.vs.immature")
merged <- merge(merged, ER.vs.immature[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".ER.vs.immature")
merged <- merge(merged, iAMP.vs.preB[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".iAMP.vs.preB")
merged <- merge(merged, DS.vs.preB[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".DS.vs.preB")
merged <- merge(merged, ER.vs.preB[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".ER.vs.preB")
merged <- merge(merged, iAMP.vs.mature[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".iAMP.vs.mature")
merged <- merge(merged, DS.vs.mature[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".DS.vs.mature")
merged <- merge(merged, ER.vs.mature[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-1):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-1):ncol(merged)], ".ER.vs.mature")

# generate unique row name for each gene
rowname <- with(merged, ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", as.character(hgnc_symbol), as.character(id)))
rowname[duplicated(rowname)] <- as.character(merged$id[duplicated(rowname)])
rownames(merged) <- rowname

cols.fc <- c("log2FoldChange.iAMP.vs.ER", "log2FoldChange.iAMP.vs.DS", "log2FoldChange.ER.vs.DS", 
		"log2FoldChange.iAMP.vs.immature", "log2FoldChange.DS.vs.immature", "log2FoldChange.ER.vs.immature", 
		"log2FoldChange.iAMP.vs.preB", "log2FoldChange.DS.vs.preB", "log2FoldChange.ER.vs.preB", 
		"log2FoldChange.iAMP.vs.mature", "log2FoldChange.DS.vs.mature", "log2FoldChange.ER.vs.mature")

cols.padj <- c("padj.iAMP.vs.ER", "padj.iAMP.vs.DS", "padj.ER.vs.DS", 
		"padj.iAMP.vs.immature", "padj.DS.vs.immature", "padj.ER.vs.immature", 
		"padj.iAMP.vs.preB", "padj.DS.vs.preB", "padj.ER.vs.preB", 
		"padj.iAMP.vs.mature", "padj.DS.vs.mature", "padj.ER.vs.mature")

# replace NAs with 0
tmp <- as.matrix(merged[,cols.fc]) ; tmp[is.na(tmp)] <- 0 ; merged[,cols.fc] <- tmp
tmp <- as.matrix(merged[,cols.padj]) ; tmp[is.na(tmp)] <- 1 ; merged[,cols.padj] <- tmp

# get top FDR for each gene
merged$padj.best <- apply(merged[,cols.padj], 1, function(x) { min(x, na.rm=T) })

# heat map function
plot.heatmap <- function(data, cexCol=0.9, cexRow=0.9, sigLevel, minFC, title="") {
	#if (!is.na(minFC)) {
	#	fc.best <- apply(data[,cols.fc], 1, function (x) { max(abs(x), na.rm=T)} )
	#	data <- data[fc.best >= minFC,]
	#}
	#data <- data[data$padj.best<sigLevel,cols.fc]
	colnames(data) <- gsub("log2FoldChange.", "", colnames(data))
	print(sprintf("Number of genes in heatmap: %d", nrow(data)))
	hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
	padjs <- as.matrix(merged[rownames(data),cols.padj])
	padjs[is.na(padjs)|padjs>sigLevel] <- NA
	padjs[!is.na(padjs)&padjs<=sigLevel] <- "*"
	heatmap.2(as.matrix(data), Colv=F, Rowv=T, dendrogram="row", trace="none", col=rev(hmcol), margin=c(10, 25), cexCol=cexCol, cexRow=cexRow, keysize=0.7, 
			  colsep=seq(1:ncol(data)), rowsep=seq(1:nrow(data)), sepcolor="grey92", sepwidth=c(0.005,0.005),
			  cellnote=padjs, notecol='white',
			  main=sprintf("%s \n(max.Padj=%.4g, min.FC=%.1f)", title, sigLevel, minFC))
	
}

# ER UP
pdf("/mnt/projects/iamp/results/genes-heatmap-ER-up.pdf", height=15, width=10)
minsig <- 1e-10 ; minfc <- 2 
plot.heatmap(merged[merged$padj.iAMP.vs.ER <= minsig & merged$padj.ER.vs.DS <= minsig & merged$log2FoldChange.iAMP.vs.ER <= -minfc & merged$log2FoldChange.ER.vs.DS >= minfc, cols.fc], cexRow=0.7, sigLevel=minsig, minFC=minfc, title="ER UP")
dev.off()

# ER DN
pdf("/mnt/projects/iamp/results/genes-heatmap-ER-dn.pdf", height=15, width=10)
minsig <- 1e-5 ; minfc <- 2 
plot.heatmap(merged[merged$padj.iAMP.vs.ER <= minsig & merged$padj.ER.vs.DS <= minsig & merged$log2FoldChange.iAMP.vs.ER >= minfc & merged$log2FoldChange.ER.vs.DS <= -minfc, cols.fc], cexRow=0.7, sigLevel=minsig, minFC=minfc, title="ER DN")
dev.off()

# iAMP UP
pdf("/mnt/projects/iamp/results/genes-heatmap-iAMP-up.pdf", height=15, width=10)
minsig <- 1e-3 ; minfc <- 1 
plot.heatmap(merged[merged$padj.iAMP.vs.ER <= minsig & merged$padj.iAMP.vs.DS <= minsig & merged$log2FoldChange.iAMP.vs.ER >= minfc & merged$log2FoldChange.iAMP.vs.DS >= minfc, cols.fc], cexRow=0.7, sigLevel=minsig, minFC=minfc, title="iAMP UP")
dev.off()

# iAMP DN
pdf("/mnt/projects/iamp/results/genes-heatmap-iAMP-dn.pdf", height=15, width=10)
minsig <- 1e-2 ; minfc <- 1 
plot.heatmap(merged[merged$padj.iAMP.vs.ER <= minsig & merged$padj.iAMP.vs.DS <= minsig & merged$log2FoldChange.iAMP.vs.ER <= -minfc & merged$log2FoldChange.iAMP.vs.DS <= -minfc, cols.fc], cexRow=0.7, sigLevel=minsig, minFC=minfc, title="iAMP DN")
dev.off()

# DS UP
pdf("/mnt/projects/iamp/results/genes-heatmap-DS-up.pdf", height=15, width=10)
minsig <- 1e-4 ; minfc <- 1 
plot.heatmap(merged[merged$padj.iAMP.vs.DS <= minsig & merged$padj.ER.vs.DS <= minsig & merged$log2FoldChange.iAMP.vs.DS <= -minfc & merged$log2FoldChange.ER.vs.DS <= -minfc, cols.fc], cexRow=0.7, sigLevel=minsig, minFC=minfc, title="DS UP")
dev.off()

# DS DN
pdf("/mnt/projects/iamp/results/genes-heatmap-DS-dn.pdf", height=15, width=10)
minsig <- 1e-3 ; minfc <- 1 
plot.heatmap(merged[merged$padj.iAMP.vs.DS <= minsig & merged$padj.ER.vs.DS <= minsig & merged$log2FoldChange.iAMP.vs.DS >= minfc & merged$log2FoldChange.ER.vs.DS >= minfc, cols.fc], cexRow=0.7, sigLevel=minsig, minFC=minfc, title="DS DN")
dev.off()
