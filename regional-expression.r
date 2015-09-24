library(Gviz)
library(GenomicRanges)
library("DESeq2")

#region <- "IKZF1"; chr <- "chr7"; start <- 50000000; end <- 50800000
#region <- "CDKN2AB"; chr <- "chr9"; start <- 21600000; end <- 22400000
#region <- "PAR1"; chr <- "chrX"; start <- 1000000; end <- 2000000
#region <- "PAX5"; chr <- "chr9"; start <- 36495626; end <- 37377382
#region <- "GABRB3"; chr <- "chr15"; start <- 25870574; end <- 27936343
#region <- "ACSM2A"; chr <- "chr16"; start <- 20318327; end <- 20643523
#region <- "BTG1"; chr <- "chr12"; start <- 91278699; end <- 92889529

# compute counts ourselves
#---
annotation <- read.delim("/mnt/projects/iamp/data/qlucore/annotations.txt")
#count.files <- list.files(path="/mnt/projects/iamp/results/htseq", pattern="*.count$", full.names=F)
#sample.names <- paste0(sub(".*C57C3ACXX_CV_([^_]+)_.*", "\\1", count.files))
#sample.table <- data.frame(name=sample.names, file=count.files)
#sample.table <- merge(sample.table, annotation[,c("Name", "Subtype")], by.x="name", by.y="Name")

#cds <- DESeqDataSetFromHTSeqCount(sampleTable = sample.table, directory="/mnt/projects/iamp/results/htseq", design=~1)
#cds <- estimateSizeFactors(cds)
#counts.norm <- as.data.frame(counts(cds, normalized=T))

# use counts from anduril pipeline
#---
counts.norm <- read.delim("/mnt/projects/iamp/results/anduril/execute/deseqExprMatrix/expr.csv", row.names=1)
counts.norm.noCD19 <- counts.norm[!names(counts.norm) %in% c("S1", "S2", "S3")] # throw out CD19 samples in this analysis

# scale data matrix
counts.scaled.noCD19 <- t(scale(t(as.matrix(counts.norm.noCD19))))
counts.scaled.withCD19 <- t(scale(t(as.matrix(counts.norm))))

# annotate genes with Ensembl biomart
#---
biomartfile <- "/mnt/projects/generic/data/ensembl/genes.GRCh37v75.biomart.RData"
load(biomartfile)

library(reshape)
genes <- genes[genes$chromosome_name %in% c(seq(1:22), "X", "Y", "MT"),]
genes$chromosome_name <- as.factor(as.character(genes$chromosome_name))
genes.dedup <- cast(genes[,c("ensembl_gene_id", "hgnc_symbol")], formula=ensembl_gene_id~., value="hgnc_symbol", fun.aggregate=function(x) { paste(x, collapse=",") })
names(genes.dedup) <- c("ensembl_gene_id", "hgnc_symbol")
genes <- merge(genes.dedup, unique(genes[,names(genes)[!names(genes) %in% "hgnc_symbol"]]), by="ensembl_gene_id")

counts.ann <- merge(genes[,c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")], counts.scaled.withCD19, by.x="ensembl_gene_id", by.y="row.names")
counts.ann$chromosome_name <- paste0("chr", counts.ann$chromosome_name)
counts.ann$hgnc_symbol <- as.character(counts.ann$hgnc_symbol)
counts.ann <- counts.ann[!duplicated(paste0(counts.ann$chromosome_name, counts.ann$start_position)),] # only one value per coordinate

gr <- makeGRangesFromDataFrame(counts.ann,
		keep.extra.columns=TRUE,
		ignore.strand=TRUE,
		seqinfo=NULL,
		seqnames.field=c("chromosome_name"),
		start.field=c("start_position"),
		end.field=c("end_position"))
gr.data <- gr[,names(mcols(gr))[!names(mcols(gr)) %in% c("ensembl_gene_id", "hgnc_symbol")]]

# mask centromeres
#{
#	centroms <- data.frame(chr=character(0), start=numeric(0), end=numeric(0), stringsAsFactors=F)
#	centroms[nrow(centroms)+1,] <- c("chr1", 121500000, 142618768)
#	centroms[nrow(centroms)+1,] <- c("chr2", 90500000, 96800000)
#	centroms[nrow(centroms)+1,] <- c("chr3", 87900000, 94225609)
#	centroms[nrow(centroms)+1,] <- c("chr4", 48200000, 52700000)
#	centroms[nrow(centroms)+1,] <- c("chr5", 46100000, 52700000)
#	centroms[nrow(centroms)+1,] <- c("chr6", 58700000, 63921350)
#	centroms[nrow(centroms)+1,] <- c("chr7", 58000000, 62574012)
#	centroms[nrow(centroms)+1,] <- c("chr8", 43100000, 48083333)
#	centroms[nrow(centroms)+1,] <- c("chr9", 47300000, 66457284)
#	centroms[nrow(centroms)+1,] <- c("chr10", 38000000, 42827313)
#	centroms[nrow(centroms)+1,] <- c("chr11", 51600000, 55700000)
#	centroms[nrow(centroms)+1,] <- c("chr12", 33300000, 38200000)
#	centroms[nrow(centroms)+1,] <- c("chr13", 16300000, 19244000)
#	centroms[nrow(centroms)+1,] <- c("chr14", 16100000, 19100000)
#	centroms[nrow(centroms)+1,] <- c("chr15", 15800000, 20178034)
#	centroms[nrow(centroms)+1,] <- c("chr16", 34600000, 46503252)
#	centroms[nrow(centroms)+1,] <- c("chr17", 22200000, 25788000)
#	centroms[nrow(centroms)+1,] <- c("chr18", 15400000, 19000000)
#	centroms[nrow(centroms)+1,] <- c("chr19", 24400000, 28600000)
#	centroms[nrow(centroms)+1,] <- c("chr20", 25752000, 29400000)
#	centroms[nrow(centroms)+1,] <- c("chr21", 10911333, 14410480)
#	centroms[nrow(centroms)+1,] <- c("chr22", 12200000, 17064000)
#	centroms[nrow(centroms)+1,] <- c("chrX", 58100000, 62983667)
#	centroms[nrow(centroms)+1,] <- c("chrY", 10037916, 14460539)
#	
#	intoNbins <- function(gr, n = 10)
#	{
#		if (any(width(gr) < n))
#			stop("all 'width(gr)' must be >= 'n'")
#		d <- width(gr) / n
#		dd <- cumsum(rep(d, each=n))
#		mask <- logical(n); mask[1] <- TRUE
#		dd <- dd - rep(dd[mask], each=n)
#		
#		starts <- round(rep(start(gr), each=n) + dd)
#		ends <- c(starts[-1], 0) - 1L
#		ends[rev(mask)] <- end(gr)
#		
#		gr <- gr[rep(seq_along(gr), each=n)]
#		ranges(gr) <- IRanges(starts, ends)
#		gr
#	}	
#
#	for (chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {
#		cen.gr <- gr.data[seqnames(gr.data)==chr,][1]
#		centrom <- centroms[centroms$chr==chr,]
#		ranges(cen.gr) <- IRanges(as.numeric(centrom$start), as.numeric(centrom$end))
#		values(cen.gr) <- as.matrix(values(cen.gr)) * 0
#		cen.gr.bins <- intoNbins(cen.gr, 200)
#		gr.data <- c(gr.data, cen.gr.bins)
#	}
#}
gr.data <- sort(gr.data)

subtypes <- as.character(annotation$Subtype[match(names(mcols(gr.data)), annotation$Name)])
subtypes[annotation$Name=="S1"] <- "immature"
subtypes[annotation$Name=="S2"] <- "preB"
subtypes[annotation$Name=="S3"] <- "mature"
subtypes.withCD19 <- factor(subtypes, levels=c("ER", "DS", "iAMP", "immature", "preB", "mature"))
subtypes.withoutCD19 <- factor(subtypes[!subtypes %in% c("immature", "preB", "mature")], levels=c("ER", "DS", "iAMP"))

# get gene models in region of interest
options(ucscChromosomeNames=FALSE)
#biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = end, name = "Genes", showId=TRUE)

# keep only longest transcript per gene
#trlen <- aggregate(width ~ gene + transcript, data=as.data.frame(biomTrack@range), FUN=sum)
#biomTrack@range <- biomTrack@range[biomTrack@range$transcript %in% trlen[ave(trlen$width, trlen$gene, FUN=max)==trlen$width,"transcript"]]
#seqlevels(ranges(biomTrack)) <- paste0("chr", seqlevels(ranges(biomTrack)))

# common region of amplification (Strefford et al., 2006)
ann.strefford <- data.frame(id=character(), seqnames=character(0), start = numeric(0), end = numeric(0), stringsAsFactors=F)
ann.strefford <- rbind(ann.strefford, data.frame(id="CRA", seqnames="chr21", start = 33192000, end = 39796000, stringsAsFactors=F)) 
ann.strefford <- rbind(ann.strefford, data.frame(id="CRD", seqnames="chr21", start = 43700000, end = 47000000, stringsAsFactors=F))

# down syndrome critical region (Olson et al., 2004); ranges from gene CBR1 to gene MX1 (Table S1)
ann.olson <- data.frame(id="DSCR", seqnames="chr21", start = 37442000, end = 42800000, stringsAsFactors=F)

# smooth with moving averages
#---
gr.ma3 <- gr.data
values(gr.ma3) <- caTools::runmean(as.matrix(values(gr.data)), 3)
names(mcols(gr.ma3)) <- names(mcols(gr.data))

gr.ma10 <- gr.data
values(gr.ma10) <- caTools::runmean(as.matrix(values(gr.data)), 10)
names(mcols(gr.ma10)) <- names(mcols(gr.data))

gr.ma50 <- gr.data
values(gr.ma50) <- caTools::runmean(as.matrix(values(gr.data)), 50)
names(mcols(gr.ma50)) <- names(mcols(gr.data))

gr.ma200 <- gr.data
values(gr.ma200) <- caTools::runmean(as.matrix(values(gr.data)), 200)
names(mcols(gr.ma200)) <- names(mcols(gr.data))

# smooth with cubic splines
# this code works, but does not really produce better results than the moving averages; maybe try loess...
# 'df' parameter controls smoothness of curve (lower is smoother)
#---
#gr.spline <- gr.data
#for (c in unique(seqnames(gr.spline))) {
#  x <- start(gr.spline[seqnames(gr.spline)==c,])
#  for (s in names(values(gr.spline))) {
#    y <- values(gr.spline[seqnames(gr.spline)==c,])[[s]]
#    smoothed <- smooth.spline(x, y, tol=1, df=10)
#    values(gr.spline[seqnames(gr.spline)==c,s]) <- smoothed$y
#  }
#}

# chromosomal plots (subtype averages)
#---
pdf("/mnt/projects/iamp/results/regional-expression.averages.pdf", width=18, height=12)
for (chr in c("chr3")) {
#for (chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {
	
	# ideogram
	itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
	
	# axis track with highlighted centromere
	cen.start <- min(start(itrack@range[itrack@range$type=="acen"]))
	cen.end <- max(end(itrack@range[itrack@range$type=="acen"]))
	gtrack <- GenomeAxisTrack(range = IRanges(start = cen.start, end = cen.end, names="CEN"), showId=TRUE)
	
	dtrack.all <- DataTrack(gr.data[,!names(mcols(gr.data)) %in% c("S1", "S2", "S3")], name="Z-score FPM", chromosome = chr, groups=subtypes.withoutCD19, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), ylim=c(-2,2))
#	dtrack.win100 <- DataTrack(gr.data, name="100 bins", chromosome = chr, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), window=100, ylim=c(-1.5,1.5))
#	dtrack.win50 <- DataTrack(gr.data, name="50 bins", chromosome = chr, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), window=50, ylim=c(-1,1))
#	dtrack.win10 <- DataTrack(gr.data, name="10 bins", chromosome = chr, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), window=10, ylim=c(-0.5,0.5), legend=TRUE)
	
	dtrack.ma3 <- DataTrack(gr.ma3[,!names(mcols(gr.ma3)) %in% c("S1", "S2", "S3")], name="mAvg (k=3)", chromosome = chr, groups=subtypes.withoutCD19, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), ylim=c(-1.5,1.5))
	dtrack.ma10 <- DataTrack(gr.ma10[,!names(mcols(gr.ma10)) %in% c("S1", "S2", "S3")], name="mAvg (k=10)", chromosome = chr, groups=subtypes.withoutCD19, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), lty=c(1,1,1), ylim=c(-1,1))
	dtrack.ma50 <- DataTrack(gr.ma50, name="mAvg (k=50)", chromosome = chr, groups=subtypes.withCD19, type = c("a", "confint", "g"), col=c("blue", "red", "orange", "black", "black", "black"), lty=c(1,1,1,1,2,3), ylim=c(-0.75,0.75))
	dtrack.ma200 <- DataTrack(gr.ma200, name="mAvg (k=200)", chromosome = chr, groups=subtypes.withCD19, type = c("a", "confint", "g"), col=c("blue", "red", "orange", "black", "black", "black"), lty=c(1,1,1,1,2,3), ylim=c(-0.5,0.5), legend=TRUE)
	#dtrack.win2 <- DataTrack(gr.data, name="2 bins", chromosome = chr, type = c("a", "confint", "g"), col=c("blue", "red", "orange"), window=2, ylim=c(-0.5,0.5), legend=TRUE)

	if (chr %in% ann.strefford$seqnames) {
		atrack.1 <- AnnotationTrack(ann.strefford, name="Strefford (2006)", col="black", cex=0.75, fontcolor.feature="black", rot.title=1, featureAnnotation = "id", chromosome = chr)
		atrack.2 <- AnnotationTrack(ann.olson, name="Olson (2004)", col="black", cex=0.75, fontcolor.feature="black", rot.title=1, featureAnnotation = "id", chromosome = chr)
		plotTracks(c(itrack, gtrack, atrack.1, atrack.2, dtrack.all, dtrack.ma3, dtrack.ma10, dtrack.ma50, dtrack.ma200), cex.title=0.7, cex.axis=0.6)	
	} else {
		plotTracks(c(itrack, gtrack, dtrack.all, dtrack.ma3, dtrack.ma10, dtrack.ma50, dtrack.ma200), cex.title=0.7, cex.axis=0.6)	
	}
}
dev.off()

# chromosomal plots (individual samples)
#---

# read cnvs
cnvs.df <- read.delim("/mnt/projects/iamp/data/iamp_cytoscan_cnvs.tsv", stringsAsFactors=F)

cnvs.df$sample[cnvs.df$sample=="1068D.cyhd.cychp"] <- "A1"
cnvs.df$sample[cnvs.df$sample=="768D.cyhd.cychp"] <- "A2"
cnvs.df$sample[cnvs.df$sample=="499D.cyhd.cychp"] <- "A3"
cnvs.df$sample[cnvs.df$sample=="1117D.cyhd.cychp"] <- "A4"
cnvs.df$sample[cnvs.df$sample=="92_D.cyhd.cychp"] <- "A5"
cnvs.df$sample[cnvs.df$sample=="1092D.cyhd.cychp"] <- "A9"
cnvs.df$sample[cnvs.df$sample=="839_R.cyhd.cychp"] <- "A10"
cnvs.df$sample[cnvs.df$sample=="456R.cyhd.cychp"] <- "A11"

cnvs.df$chr <- paste0("chr", cnvs.df$chr)
cnvs.df$integerCN <- cnvs.df$cn
cnvs.df$integerCN[cnvs.df$integerCN>2 & cnvs.df$integerCN<3] <- 3
cnvs.df$integerCN[cnvs.df$integerCN<2 & cnvs.df$integerCN>1] <- 1
cnvs.df$integerCN <- round(cnvs.df$integerCN)

cnvs.gr <- as(cnvs.df, "GRanges")

# find common region of amplification (CRA)
cnvs.gr.cra <- NA
cnvs.samples <- unique(cnvs.gr$sample)
for (i in 1:length(cnvs.samples)) {
  iamps <- cnvs.gr[cnvs.gr$sample == cnvs.samples[i] & seqnames(cnvs.gr)=="chr21" & cnvs.gr$integerCN > 2,]
  if (i == 1) {
    cnvs.gr.cra <- iamps
  } else {
    cnvs.gr.cra <- intersect(cnvs.gr.cra, iamps)
  }
}

counts.ann.withCD19 <- merge(genes[,c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")], counts.scaled.withCD19, by.x="ensembl_gene_id", by.y="row.names")
counts.ann.withCD19$chromosome_name <- paste0("chr", counts.ann.withCD19$chromosome_name)

gr.withCD19 <- makeGRangesFromDataFrame(counts.ann.withCD19,
		keep.extra.columns=TRUE,
		ignore.strand=TRUE,
		seqinfo=NULL,
		seqnames.field=c("chromosome_name"),
		start.field=c("start_position"),
		end.field=c("end_position"))
gr.withCD19 <- sort(gr.withCD19)
gr.withCD19 <- gr.withCD19[,names(mcols(gr.withCD19))[!names(mcols(gr.withCD19)) %in% c("ensembl_gene_id", "hgnc_symbol")]]

gr.withCD19.ma200 <- gr.withCD19
values(gr.withCD19.ma200) <- caTools::runmean(as.matrix(values(gr.withCD19)), 200)
names(mcols(gr.withCD19.ma200)) <- names(mcols(gr.withCD19))

gr.iamp <- gr.withCD19.ma200[,grep("^A", names(mcols(gr.withCD19.ma200)), value=T)]
gr.er <- gr.withCD19.ma200[,grep("^C", names(mcols(gr.withCD19.ma200)), value=T)]
gr.pc <- gr.withCD19.ma200[,grep("^D", names(mcols(gr.withCD19.ma200)), value=T)]
gr.cd19 <- gr.withCD19.ma200[,grep("^S", names(mcols(gr.withCD19.ma200)), value=T)]

# DEGs
degs <- read.delim("/mnt/projects/iamp/results/anduril/execute/degTable_iAMPvsNoniAMP/table.csv", stringsAsFactors = F)
degs$chr <- paste0("chr", degs$chr)
degs$status <- ifelse(degs$status > 0, "UP", "DOWN")
gr.degs <- makeGRangesFromDataFrame(degs, keep.extra.columns=TRUE, ignore.strand=TRUE, seqinfo=NULL)

pdf("/mnt/projects/iamp/results/regional-expression.samples.pdf", width=18, height=12)
#for (chr in c("chr21")) {
for (chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {
	itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
	
	# axis track with highlighted centromere
	cen.start <- min(start(itrack@range[itrack@range$type=="acen"]))
	cen.end <- max(end(itrack@range[itrack@range$type=="acen"]))
	gtrack <- GenomeAxisTrack(range = IRanges(start = cen.start, end = cen.end, names="CEN"), showId=TRUE)
	
	cntrack <- AnnotationTrack(cnvs.gr, chromosome=chr, name="CNA", rot.title=1, id=cnvs.gr$sample, shape = "box", groupAnnotation = "id", stackHeight=0.9, cex.group=0.8, fontcolor.feature="black", group=cnvs.gr$sample, min.width=1, col.line = "lightgray", col = NULL)
	feature(cntrack) <- paste0("CN", cnvs.gr$integerCN)
	
	dtrack.er <- DataTrack(gr.er, name="ER", chromosome = chr, type = c("a", "g"), groups=names(mcols(gr.er)), ylim=c(-0.75,0.75), legend=TRUE)
	dtrack.pc <- DataTrack(gr.pc, name="DS", chromosome = chr, type = c("a", "g"), groups=names(mcols(gr.pc)), ylim=c(-0.75,0.75), legend=TRUE)
	dtrack.iamp <- DataTrack(gr.iamp, name="iAMP21", chromosome = chr, type = c("a", "g"), groups=names(mcols(gr.iamp)), ylim=c(-0.75,0.75), legend=TRUE)
	dtrack.cd19 <- DataTrack(gr.cd19, name="CD19", chromosome = chr, type = c("a", "g"), groups=names(mcols(gr.cd19)), ylim=c(-0.75,0.75), legend=TRUE)

	atrack.degs <- AnnotationTrack(gr.degs, chromosome=chr, name="DEGs", id=gr.degs$Gene, feature=gr.degs$status, rot.title=1, shape = "box", group=gr.degs$Gene, groupAnnotation = "id", cex.group=0.7, stackHeight=1, fontcolor.feature="black", min.width=1, col = NULL)
	atrack.genes <- AnnotationTrack(gr, chromosome=chr, name="Genes", id=gr$Ensembl, rot.title=1, shape = "box", stacking = "dense", min.width=1, fill = "lightgray", col = "black", alpha = 0.7, alpha.title = 1)
	
	if (chr == "chr21") {
	  atrack.cra <- AnnotationTrack(cnvs.gr.cra, chromosome=chr, name="CRA", id="Shared", feature="CN4", shape = "box", groupAnnotation = "id", cex.group=0.8, fontcolor.feature="black", group="CRA", min.width=1, col.line = "lightgray", col = NULL, rot.title=1)
	  atrack.strefford <- AnnotationTrack(ann.strefford, name="Streff.", col="black", cex=0.75, fontcolor.feature="black", rot.title=1, featureAnnotation = "id", chromosome = chr)
	  sizes <- c(0.03, 0.02, 0.03, 0.14, 0.02, 0.02, 0.1, 0.16, 0.16, 0.16, 0.16)
	  plotTracks(c(itrack, gtrack, atrack.strefford, cntrack, atrack.cra, atrack.genes, atrack.degs, dtrack.iamp, dtrack.pc, dtrack.er, dtrack.cd19), cex.title=0.9, cex.axis=0.7, sizes=sizes, CN0="darkblue", CN1="lightblue", CN2="white", CN3="red", CN4="darkred", CRA="darkred", DOWN="blue", UP="red", collapse=FALSE)	
	} else {
		sizes <- c(0.03, 0.05, 0.07, 0.03, 0.08, 0.185, 0.185, 0.185, 0.185)
		plotTracks(c(itrack, gtrack, cntrack, atrack.genes, atrack.degs, dtrack.iamp, dtrack.pc, dtrack.er, dtrack.cd19), cex.title=0.9, cex.axis=0.7, sizes=sizes, CN0="darkblue", CN1="lightblue", CN2="white", CN3="red", CN4="darkred", collapse=FALSE)	
	}
}
dev.off()

# boxplots of individual genes
#---
library(lattice)

counts.melt <- melt(as.matrix(counts.norm))
names(counts.melt) <- c("ensembl", "sample", "count")
counts.melt <- merge(counts.melt, annotation[,c("Name", "Subtype")], by.x="sample", by.y="Name")
counts.melt$Subtype <- factor(as.character(counts.melt$Subtype), levels=c("CD19", "ER", "DS", "iAMP"))
counts.melt <- merge(genes[,c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")], counts.melt, by.x="ensembl_gene_id", by.y="ensembl")
counts.melt$title <- ifelse(counts.melt$hgnc_symbol != "", as.character(counts.melt$hgnc_symbol), counts.melt$ensembl_gene_id)
counts.melt$title <- factor(as.character(counts.melt$title), levels=unique(counts.melt$title[order(counts.melt$start)]))

# CRA
rgenes <- counts.melt[counts.melt$chromosome_name=="21" & counts.melt$start_position >= 33192000 & counts.melt$end_position <= 39796000,]
pdf("/mnt/projects/iamp/results/regional-expression-CRA.pdf", width=14, height=10)
print(bwplot(count~Subtype | title, data=rgenes, 
				layout=c(17, (length(unique(rgenes$title))+16) %/% 17),
				par.strip.text=list(cex=0.5),
				main="Common region of amplification (Strefford et al., 2006)\nchr21:33,192,000-39,796,000",
				notch=FALSE,
				as.table=TRUE,
				ylab="log(FPM)",
				scales=list(x=list(rot=90)),
				par.settings = list(box.umbrella=list(col="black"), box.rectangle = list(col="black")), 
				panel=function(x,y,...){
					panel.grid()
					panel.bwplot(x,y,pch="|",do.out=FALSE, ...)
					panel.stripplot(x,y,jitter.data=TRUE,factor=0.8,pch=19,cex=0.3,...)
				}))
dev.off()

# CRD
rgenes <- counts.melt[counts.melt$chromosome_name=="21" & counts.melt$start_position >= 43700000 & counts.melt$end_position <= 47000000,]
pdf("/mnt/projects/iamp/results/regional-expression-CRD.pdf", width=14, height=10)
print(bwplot(count~Subtype | title, data=rgenes, 
				layout=c(17, (length(unique(rgenes$title))+16) %/% 17),
				par.strip.text=list(cex=0.5),
				main="Common region of deletion (Strefford et al., 2006)\nchr21:43,700,000-47,000,000",
				notch=FALSE,
				as.table=TRUE,
				ylab="log(FPM)",
				scales=list(x=list(rot=90)),
				par.settings = list(box.umbrella=list(col="black"), box.rectangle = list(col="black")), 
				panel=function(x,y,...){
					panel.grid()
					panel.bwplot(x,y,pch="|",do.out=FALSE, ...)
					panel.stripplot(x,y,jitter.data=TRUE,factor=0.8,pch=19,cex=0.3,...)
				}))
dev.off()

# DSCR
rgenes <- counts.melt[counts.melt$chromosome_name=="21" & counts.melt$start_position >= 37442000 & counts.melt$end_position <= 42800000,]
pdf("/mnt/projects/iamp/results/regional-expression-DSCR.pdf", width=14, height=10)
print(bwplot(count~Subtype | title, data=rgenes, 
				layout=c(17, (length(unique(rgenes$title))+16) %/% 17),
				par.strip.text=list(cex=0.5),
				main="Down syndrome critical region (Olson et al., 2004)\nchr21:37,442,000-42,800,000",
				notch=FALSE,
				as.table=TRUE,
				ylab="log(FPM)",
				scales=list(x=list(rot=90)),
				par.settings = list(box.umbrella=list(col="black"), box.rectangle = list(col="black")), 
				panel=function(x,y,...){
					panel.grid()
					panel.bwplot(x,y,pch="|",do.out=FALSE, ...)
					panel.stripplot(x,y,jitter.data=TRUE,factor=0.8,pch=19,cex=0.3,...)
				}))
dev.off()

stop("OK")

# plot genes in region plots
chr <- "chr21"; start <- 33192000; end <- 34000000
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
gtrack <- GenomeAxisTrack()

# gene track with only longest isoforms
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = end, name = "Genes", showId=TRUE)
trlen <- aggregate(width ~ gene + transcript, data=as.data.frame(biomTrack@range), FUN=sum)
biomTrack@range <- biomTrack@range[biomTrack@range$transcript %in% trlen[ave(trlen$width, trlen$gene, FUN=max)==trlen$width,"transcript"] & biomTrack@range$gene %in% rownames(counts.norm)]

plotTracks(c(itrack, gtrack, biomTrack), cex.title=0.9)	

#plotTracks(c(list(itrack, gtrack, biomTrack), cntracks), from=start, to=end, title.width=1.8, CN0="darkred", "CN1"="red", "CN2"="white", "CN3"="lightblue", "CN4"="darkblue")

#cntrack@range@elementMetadata@listData$id <- gr$sample.name

#pdf("/mnt/projects/p2ry8-crlf2/results/exomeCopy/IKZF1.pdf")
#plotTracks(list(itrack, gtrack, biomTrack, cntrack), from=start, to=end)
#dev.off()
