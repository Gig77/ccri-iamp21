options(warn=1)

iAMP.vs.PC <- read.delim("~/iamp/results/deseq/iAMP-vs-PC.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.PC <- iAMP.vs.PC[!is.na(iAMP.vs.PC$pvalue),]
iAMP.vs.PC <- iAMP.vs.PC[ave(iAMP.vs.PC$pvalue, iAMP.vs.PC$id, FUN=min) == iAMP.vs.PC$pvalue,]

iAMP.vs.ER <- read.delim("~/iamp/results/deseq/iAMP-vs-ER.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.ER <- iAMP.vs.ER[!is.na(iAMP.vs.ER$pvalue),]
iAMP.vs.ER <- iAMP.vs.ER[ave(iAMP.vs.ER$pvalue, iAMP.vs.ER$id, FUN=min) == iAMP.vs.ER$pvalue,]

ER.vs.PC <- read.delim("~/iamp/results/deseq/ER-vs-PC.tsv", check.names=F, stringsAsFactors=F)
ER.vs.PC <- ER.vs.PC[!is.na(ER.vs.PC$pvalue),]
ER.vs.PC <- ER.vs.PC[ave(ER.vs.PC$pvalue, ER.vs.PC$id, FUN=min) == ER.vs.PC$pvalue,]

iAMP.vs.immature <- read.delim("~/iamp/results/deseq/iAMP-vs-immature.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.immature <- iAMP.vs.immature[!is.na(iAMP.vs.immature$pvalue),]
iAMP.vs.immature <- iAMP.vs.immature[ave(iAMP.vs.immature$pvalue, iAMP.vs.immature$id, FUN=min) == iAMP.vs.immature$pvalue,]

iAMP.vs.preB <- read.delim("~/iamp/results/deseq/iAMP-vs-preB.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.preB <- iAMP.vs.preB[!is.na(iAMP.vs.preB$pvalue),]
iAMP.vs.preB <- iAMP.vs.preB[ave(iAMP.vs.preB$pvalue, iAMP.vs.preB$id, FUN=min) == iAMP.vs.preB$pvalue,]

iAMP.vs.mature <- read.delim("~/iamp/results/deseq/iAMP-vs-mature.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.mature <- iAMP.vs.mature[!is.na(iAMP.vs.mature$pvalue),]
iAMP.vs.mature <- iAMP.vs.mature[ave(iAMP.vs.mature$pvalue, iAMP.vs.mature$id, FUN=min) == iAMP.vs.mature$pvalue,]

PC.vs.immature <- read.delim("~/iamp/results/deseq/PC-vs-immature.tsv", check.names=F, stringsAsFactors=F)
PC.vs.immature <- PC.vs.immature[!is.na(PC.vs.immature$pvalue),]
PC.vs.immature <- PC.vs.immature[ave(PC.vs.immature$pvalue, PC.vs.immature$id, FUN=min) == PC.vs.immature$pvalue,]

PC.vs.preB <- read.delim("~/iamp/results/deseq/PC-vs-preB.tsv", check.names=F, stringsAsFactors=F)
PC.vs.preB <- PC.vs.preB[!is.na(PC.vs.preB$pvalue),]
PC.vs.preB <- PC.vs.preB[ave(PC.vs.preB$pvalue, PC.vs.preB$id, FUN=min) == PC.vs.preB$pvalue,]

PC.vs.mature <- read.delim("~/iamp/results/deseq/PC-vs-mature.tsv", check.names=F, stringsAsFactors=F)
PC.vs.mature <- PC.vs.mature[!is.na(PC.vs.mature$pvalue),]
PC.vs.mature <- PC.vs.mature[ave(PC.vs.mature$pvalue, PC.vs.mature$id, FUN=min) == PC.vs.mature$pvalue,]

ER.vs.immature <- read.delim("~/iamp/results/deseq/ER-vs-immature.tsv", check.names=F, stringsAsFactors=F)
ER.vs.immature <- ER.vs.immature[!is.na(ER.vs.immature$pvalue),]
ER.vs.immature <- ER.vs.immature[ave(ER.vs.immature$pvalue, ER.vs.immature$id, FUN=min) == ER.vs.immature$pvalue,]

ER.vs.preB <- read.delim("~/iamp/results/deseq/ER-vs-preB.tsv", check.names=F, stringsAsFactors=F)
ER.vs.preB <- ER.vs.preB[!is.na(ER.vs.preB$pvalue),]
ER.vs.preB <- ER.vs.preB[ave(ER.vs.preB$pvalue, ER.vs.preB$id, FUN=min) == ER.vs.preB$pvalue,]

ER.vs.mature <- read.delim("~/iamp/results/deseq/ER-vs-mature.tsv", check.names=F, stringsAsFactors=F)
ER.vs.mature <- ER.vs.mature[!is.na(ER.vs.mature$pvalue),]
ER.vs.mature <- ER.vs.mature[ave(ER.vs.mature$pvalue, ER.vs.mature$id, FUN=min) == ER.vs.mature$pvalue,]

# merge the three data sets, add suffixes
merge.cols <- c("id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "description")
value.cols <- c(1:6,8:10,13,14)

merged <- merge(iAMP.vs.PC[,value.cols], iAMP.vs.ER[,value.cols], by=merge.cols, suffixes=c(".iAMP.vs.PC", ".iAMP.vs.ER"), all=T)
merged <- merge(merged, ER.vs.PC[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.PC")

merged <- merge(merged, iAMP.vs.immature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".iAMP.vs.immature")
merged <- merge(merged, iAMP.vs.preB[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".iAMP.vs.preB")
merged <- merge(merged, iAMP.vs.mature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".iAMP.vs.mature")

merged <- merge(merged, PC.vs.immature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".PC.vs.immature")
merged <- merge(merged, PC.vs.preB[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".PC.vs.preB")
merged <- merge(merged, PC.vs.mature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".PC.vs.mature")

merged <- merge(merged, ER.vs.immature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.immature")
merged <- merge(merged, ER.vs.preB[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.preB")
merged <- merge(merged, ER.vs.mature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.mature")

# merge miR annotation: validated targets
library(reshape)
validated <- read.delim("~/iamp/results/miRecords.v4.validatedTargets.txt")
validated.collapsed <- cast(validated, formula = symbol ~ ., value="miRNA_mature_ID", fun.aggregate=function(x) { paste(unique(x), collapse=";") }, fill="")
names(validated.collapsed)[2] <- "miR.validated"
merged <- merge(merged, validated.collapsed, by.x="hgnc_symbol", by.y="symbol", all.x=T)

# merge miR annotation: predicted targets
predicted <- read.delim("~/iamp/results/miRecords.v4.predictedTargets.by5programs.txt")
predicted.collapsed <- cast(predicted, formula = Symbol ~ ., value="miRNA.ID", fun.aggregate=function(x) { paste(unique(x), collapse=";") }, fill="")
names(predicted.collapsed)[2] <- "miR.predicted"
merged <- merge(merged, predicted.collapsed, by.x="hgnc_symbol", by.y="Symbol", all.x=T)

# reorder columns
merged.reordered <- merged[,c(1:6,67,68,7:66)]

# write table
write.table(merged.reordered, file="~/iamp/results/deseq/iAMP21.diff-exp-genes.deseq2.merged.pairwise.tsv.part", row.names=F, col.names=T, quote=F, sep="\t")
