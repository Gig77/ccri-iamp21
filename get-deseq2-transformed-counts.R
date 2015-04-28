options(warn=1)
library("DESeq2")

# build dataframe with count files
files <- list.files(path="~/iamp/results/htseq/", pattern=".count$")
names <- sub(".*CV_([^_]+).*", "\\1", files)
samples <- data.frame(name=names, file=files, stringsAsFactors=F)

# transform counts into normalized values
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="~/iamp/results/htseq", design=~1)

# regularized log transformation
rld <- rlog(cds)
rlogMat <- assay(rld)

# variance stabilizing transformation
#vsd <- varianceStabilizingTransformation(cds)
#vstMat <- assay(vsd)

# annotate gene names
#rlogMat <- read.delim("~/fikret/results/deseq/normalized-counts.deseq2.rlog.tsv", check.names=F)
load("~/generic/data/ensembl/genes.GRCh37v75.biomart.RData")
rlogMat <- as.data.frame(rlogMat)
rlogMat$ensembl_gene_id <- rownames(rlogMat)
rlogMat.hgnc <- merge(rlogMat, genes[,1:8], all.x=T)

# write table
n <- names(rlogMat.hgnc)
write.table(rlogMat.hgnc[n[c(1,(length(n)-6):length(n),2:(length(n)-7))]], file="~/iamp/results/qlucore/normalized-counts.deseq2.rlog.tsv", col.names=T, row.names=F, sep="\t", quote=F)

