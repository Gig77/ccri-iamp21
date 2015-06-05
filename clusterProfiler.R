# install
source("http://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
biocLite("org.Hs.eg.db")

library(clusterProfiler)

# example data
library(DOSE)
data(geneList)

# get genes IDs and their log2 fold changes
deseq <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-DS.tsv")

# merge Entrez gene ID
library(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
bm <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), mart=mart)
deseq <- unique(merge(deseq, bm, by.x="id", by.y="ensembl_gene_id", all.x = TRUE))

# GO enrichment
library(GO.db)

genes.diff <-unique(deseq$entrezgene[!is.na(deseq$entrezgene) & abs(deseq$log2FoldChange) >= 2]) # differentially expressed genes
genes.bg <- unique(deseq$entrezgene[!is.na(deseq$entrezgene) & deseq$baseMean>0]) # all expresed genes

ego <- enrichGO(
  gene = genes.diff,
  universe = genes.bg,
  organism = "human",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(summary(ego))
                
# GO gene set enrichment
valid <- !is.na(deseq$pvalue) & !is.na(deseq$log2FoldChange) & !is.na(deseq$entrezgene) & !duplicated(deseq$entrezgene)
genes <- deseq$log2FoldChange[valid]
names(genes) <- deseq$entrezgene[valid]
genes <- sort(genes, decreasing=TRUE)

# NOTE: clusterProfiler produces non-sensical results, because all p-values are the same (0.009999)?
# Also, if I inject 200 genes with identical BP GO term GO:0007507 with the two lines below, this gene set does not come out as significant
#gogenes <- getBM(attributes=c('entrezgene', 'go_id'), filters = 'go_id', values = 'GO:0007507', mart = mart)
#names(genes[1:200]) <- gogenes

ego2 <- gseGO(
  geneList = genes,
  organism = "human",
  ont = "BP",
  nPerm = 100,
  minGSSize = 10,
  pvalueCutoff = 0.1,
  verbose = TRUE)
result <- summary(ego2)
head(result[order(result$pvalue),], n=20)

