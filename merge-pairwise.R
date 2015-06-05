options(warn=1)

iAMP.vs.DS <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-DS.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.DS <- iAMP.vs.DS[!is.na(iAMP.vs.DS$pvalue),]
iAMP.vs.DS <- iAMP.vs.DS[ave(iAMP.vs.DS$pvalue, iAMP.vs.DS$id, FUN=min) == iAMP.vs.DS$pvalue,]

iAMP.vs.ER <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-ER.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.ER <- iAMP.vs.ER[!is.na(iAMP.vs.ER$pvalue),]
iAMP.vs.ER <- iAMP.vs.ER[ave(iAMP.vs.ER$pvalue, iAMP.vs.ER$id, FUN=min) == iAMP.vs.ER$pvalue,]

ER.vs.DS <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-DS.tsv", check.names=F, stringsAsFactors=F)
ER.vs.DS <- ER.vs.DS[!is.na(ER.vs.DS$pvalue),]
ER.vs.DS <- ER.vs.DS[ave(ER.vs.DS$pvalue, ER.vs.DS$id, FUN=min) == ER.vs.DS$pvalue,]

iAMP.vs.immature <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-immature.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.immature <- iAMP.vs.immature[!is.na(iAMP.vs.immature$pvalue),]
iAMP.vs.immature <- iAMP.vs.immature[ave(iAMP.vs.immature$pvalue, iAMP.vs.immature$id, FUN=min) == iAMP.vs.immature$pvalue,]

iAMP.vs.preB <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-preB.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.preB <- iAMP.vs.preB[!is.na(iAMP.vs.preB$pvalue),]
iAMP.vs.preB <- iAMP.vs.preB[ave(iAMP.vs.preB$pvalue, iAMP.vs.preB$id, FUN=min) == iAMP.vs.preB$pvalue,]

iAMP.vs.mature <- read.delim("/mnt/projects/iamp/results/deseq/iAMP-vs-mature.tsv", check.names=F, stringsAsFactors=F)
iAMP.vs.mature <- iAMP.vs.mature[!is.na(iAMP.vs.mature$pvalue),]
iAMP.vs.mature <- iAMP.vs.mature[ave(iAMP.vs.mature$pvalue, iAMP.vs.mature$id, FUN=min) == iAMP.vs.mature$pvalue,]

DS.vs.immature <- read.delim("/mnt/projects/iamp/results/deseq/DS-vs-immature.tsv", check.names=F, stringsAsFactors=F)
DS.vs.immature <- DS.vs.immature[!is.na(DS.vs.immature$pvalue),]
DS.vs.immature <- DS.vs.immature[ave(DS.vs.immature$pvalue, DS.vs.immature$id, FUN=min) == DS.vs.immature$pvalue,]

DS.vs.preB <- read.delim("/mnt/projects/iamp/results/deseq/DS-vs-preB.tsv", check.names=F, stringsAsFactors=F)
DS.vs.preB <- DS.vs.preB[!is.na(DS.vs.preB$pvalue),]
DS.vs.preB <- DS.vs.preB[ave(DS.vs.preB$pvalue, DS.vs.preB$id, FUN=min) == DS.vs.preB$pvalue,]

DS.vs.mature <- read.delim("/mnt/projects/iamp/results/deseq/DS-vs-mature.tsv", check.names=F, stringsAsFactors=F)
DS.vs.mature <- DS.vs.mature[!is.na(DS.vs.mature$pvalue),]
DS.vs.mature <- DS.vs.mature[ave(DS.vs.mature$pvalue, DS.vs.mature$id, FUN=min) == DS.vs.mature$pvalue,]

ER.vs.immature <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-immature.tsv", check.names=F, stringsAsFactors=F)
ER.vs.immature <- ER.vs.immature[!is.na(ER.vs.immature$pvalue),]
ER.vs.immature <- ER.vs.immature[ave(ER.vs.immature$pvalue, ER.vs.immature$id, FUN=min) == ER.vs.immature$pvalue,]

ER.vs.preB <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-preB.tsv", check.names=F, stringsAsFactors=F)
ER.vs.preB <- ER.vs.preB[!is.na(ER.vs.preB$pvalue),]
ER.vs.preB <- ER.vs.preB[ave(ER.vs.preB$pvalue, ER.vs.preB$id, FUN=min) == ER.vs.preB$pvalue,]

ER.vs.mature <- read.delim("/mnt/projects/iamp/results/deseq/ER-vs-mature.tsv", check.names=F, stringsAsFactors=F)
ER.vs.mature <- ER.vs.mature[!is.na(ER.vs.mature$pvalue),]
ER.vs.mature <- ER.vs.mature[ave(ER.vs.mature$pvalue, ER.vs.mature$id, FUN=min) == ER.vs.mature$pvalue,]

# merge the three data sets, add suffixes
merge.cols <- c("id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "description")
value.cols <- c(1:6,8:10,13,14)

merged <- merge(iAMP.vs.DS[,value.cols], iAMP.vs.ER[,value.cols], by=merge.cols, suffixes=c(".iAMP.vs.DS", ".iAMP.vs.ER"), all=T)
merged <- merge(merged, ER.vs.DS[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.DS")

merged <- merge(merged, iAMP.vs.immature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".iAMP.vs.immature")
merged <- merge(merged, iAMP.vs.preB[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".iAMP.vs.preB")
merged <- merge(merged, iAMP.vs.mature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".iAMP.vs.mature")

merged <- merge(merged, DS.vs.immature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".DS.vs.immature")
merged <- merge(merged, DS.vs.preB[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".DS.vs.preB")
merged <- merge(merged, DS.vs.mature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".DS.vs.mature")

merged <- merge(merged, ER.vs.immature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.immature")
merged <- merge(merged, ER.vs.preB[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.preB")
merged <- merge(merged, ER.vs.mature[,value.cols], by=merge.cols, all=T)
names(merged)[(ncol(merged)-4):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-4):ncol(merged)], ".ER.vs.mature")

# merge miR annotation: validated targets
library(reshape)
validated <- read.delim("/mnt/projects/iamp/data/miRecords/miRecords.validatedTargets.txt", stringsAsFactors = FALSE)
validated.collapsed <- cast(validated, formula = symbol ~ ., value="miRNA_mature_ID", fun.aggregate=function(x) { paste(unique(x), collapse=";") }, fill="")
names(validated.collapsed)[2] <- "miR.validated"
merged <- merge(merged, validated.collapsed, by.x="hgnc_symbol", by.y="symbol", all.x=T)

# merge miR annotation: predicted targets
predicted <- read.delim("/mnt/projects/iamp/data/miRecords/miRecords.predictedTargets.txt", stringsAsFactors = FALSE)
predicted <- predicted[predicted$databases >= 5,]
predicted.collapsed <- cast(predicted, formula = Symbol ~ ., value="miRNA.ID", fun.aggregate=function(x) { paste(unique(x), collapse=";") }, fill="")
names(predicted.collapsed)[2] <- "miR.predicted"
merged <- merge(merged, predicted.collapsed, by.x="hgnc_symbol", by.y="Symbol", all.x=T)

# reorder columns
merged.reordered <- merged[,c(1:6,67,68,7:66)]

# write table
write.table(merged.reordered, file="/mnt/projects/iamp/results/deseq/iAMP21.diff-exp-genes.deseq2.merged.pairwise.tsv.part", row.names=F, col.names=T, quote=F, sep="\t")
