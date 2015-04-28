options(warn=1)
library(VennDiagram)

iAMP.vs.PC <- read.delim("~/iamp/results/deseq/iAMP-vs-PC.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.ER <- read.delim("~/iamp/results/deseq/iAMP-vs-ER.tsv", check.names=F, stringsAsFactors=F)
ER.vs.PC <- read.delim("~/iamp/results/deseq/ER-vs-PC.tsv", check.names=F, stringsAsFactors=F)

set1 <- with(iAMP.vs.PC, iAMP.vs.PC[!is.na(padj) & padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1), "id"])
set2 <- with(iAMP.vs.ER, iAMP.vs.ER[!is.na(padj) & padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1), "id"])
set3 <- with(ER.vs.PC, ER.vs.PC[!is.na(padj) & padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1), "id"])

pdf("~/iamp/results/venn.pdf")
grid.draw(draw.triple.venn(length(set1), length(set2), length(set3), length(intersect(set1, set2)), length(intersect(set2, set3)), length(intersect(set1, set3)), length(intersect(intersect(set1, set2), set3)), category=c(sprintf("iAMP vs. PC\n(%d)", length(set1)), sprintf("iAMP vs. ER\n(%d)", length(set2)), sprintf("ER vs. PC (%d)", length(set3)))))
dev.off()
