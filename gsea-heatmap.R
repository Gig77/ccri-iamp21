library(RColorBrewer)
library("gplots")

iAMP.vs.ER.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_ER.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_ER.gsea", pattern="gsea_report_for_na_pos.*xls")))
iAMP.vs.ER.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_ER.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_ER.gsea", pattern="gsea_report_for_na_neg.*xls")))
iAMP.vs.PC.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_PC.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_PC.gsea", pattern="gsea_report_for_na_pos.*xls")))
iAMP.vs.PC.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_PC.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_PC.gsea", pattern="gsea_report_for_na_neg.*xls")))
ER.vs.PC.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_PC.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_PC.gsea", pattern="gsea_report_for_na_pos.*xls")))
ER.vs.PC.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_PC.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_PC.gsea", pattern="gsea_report_for_na_neg.*xls")))
iAMP.vs.immature.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_immature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_immature.gsea", pattern="gsea_report_for_na_pos.*xls")))
iAMP.vs.immature.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_immature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_immature.gsea", pattern="gsea_report_for_na_neg.*xls")))
PC.vs.immature.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/PC_vs_immature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/PC_vs_immature.gsea", pattern="gsea_report_for_na_pos.*xls")))
PC.vs.immature.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/PC_vs_immature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/PC_vs_immature.gsea", pattern="gsea_report_for_na_neg.*xls")))
ER.vs.immature.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_immature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_immature.gsea", pattern="gsea_report_for_na_pos.*xls")))
ER.vs.immature.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_immature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_immature.gsea", pattern="gsea_report_for_na_neg.*xls")))
iAMP.vs.preB.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_preB.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_preB.gsea", pattern="gsea_report_for_na_pos.*xls")))
iAMP.vs.preB.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_preB.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_preB.gsea", pattern="gsea_report_for_na_neg.*xls")))
PC.vs.preB.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/PC_vs_preB.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/PC_vs_preB.gsea", pattern="gsea_report_for_na_pos.*xls")))
PC.vs.preB.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/PC_vs_preB.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/PC_vs_preB.gsea", pattern="gsea_report_for_na_neg.*xls")))
ER.vs.preB.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_preB.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_preB.gsea", pattern="gsea_report_for_na_pos.*xls")))
ER.vs.preB.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_preB.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_preB.gsea", pattern="gsea_report_for_na_neg.*xls")))
iAMP.vs.mature.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_mature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_mature.gsea", pattern="gsea_report_for_na_pos.*xls")))
iAMP.vs.mature.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/iAMP_vs_mature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/iAMP_vs_mature.gsea", pattern="gsea_report_for_na_neg.*xls")))
PC.vs.mature.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/PC_vs_mature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/PC_vs_mature.gsea", pattern="gsea_report_for_na_pos.*xls")))
PC.vs.mature.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/PC_vs_mature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/PC_vs_mature.gsea", pattern="gsea_report_for_na_neg.*xls")))
ER.vs.mature.up <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_mature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_mature.gsea", pattern="gsea_report_for_na_pos.*xls")))
ER.vs.mature.down <- read.delim(paste0("/mnt/projects/iamp/results/gsea/ER_vs_mature.gsea/", list.files(path="/mnt/projects/iamp/results/gsea/ER_vs_mature.gsea", pattern="gsea_report_for_na_neg.*xls")))

cols.keep <- c(1, 5:11)
cols.merge <- c(1)
merged <- merge(iAMP.vs.ER.up[,cols.keep], iAMP.vs.ER.down[,cols.keep], by=cols.merge, suffixes=c(".iAMP.vs.ER.up", ".iAMP.vs.ER.down"), all=T)
merged <- merge(merged, iAMP.vs.PC.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.PC.up")
merged <- merge(merged, iAMP.vs.PC.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.PC.down")
merged <- merge(merged, ER.vs.PC.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.PC.up")
merged <- merge(merged, ER.vs.PC.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.PC.down")
merged <- merge(merged, iAMP.vs.immature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.immature.up")
merged <- merge(merged, iAMP.vs.immature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.immature.down")
merged <- merge(merged, PC.vs.immature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".PC.vs.immature.up")
merged <- merge(merged, PC.vs.immature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".PC.vs.immature.down")
merged <- merge(merged, ER.vs.immature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.immature.up")
merged <- merge(merged, ER.vs.immature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.immature.down")
merged <- merge(merged, iAMP.vs.preB.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.preB.up")
merged <- merge(merged, iAMP.vs.preB.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.preB.down")
merged <- merge(merged, PC.vs.preB.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".PC.vs.preB.up")
merged <- merge(merged, PC.vs.preB.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".PC.vs.preB.down")
merged <- merge(merged, ER.vs.preB.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.preB.up")
merged <- merge(merged, ER.vs.preB.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.preB.down")
merged <- merge(merged, iAMP.vs.mature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.mature.up")
merged <- merge(merged, iAMP.vs.mature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.mature.down")
merged <- merge(merged, PC.vs.mature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".PC.vs.mature.up")
merged <- merge(merged, PC.vs.mature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".PC.vs.mature.down")
merged <- merge(merged, ER.vs.mature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.mature.up")
merged <- merge(merged, ER.vs.mature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.mature.down")
rownames(merged) <- merged$NAME

# get top FDR for gene set
merged$FDR.best <- mapply(pmin, merged$FDR.q.val.iAMP.vs.ER.up, merged$FDR.q.val.iAMP.vs.ER.down, merged$FDR.q.val.iAMP.vs.PC.up, merged$FDR.q.val.iAMP.vs.PC.down, merged$FDR.q.val.ER.vs.PC.up, merged$FDR.q.val.ER.vs.PC.down, 
		                        merged$FDR.q.val.iAMP.vs.immature.up, merged$FDR.q.val.iAMP.vs.immature.down, merged$FDR.q.val.PC.vs.immature.up, merged$FDR.q.val.PC.vs.immature.down, merged$FDR.q.val.ER.vs.immature.up, merged$FDR.q.val.ER.vs.immature.down, 
								merged$FDR.q.val.iAMP.vs.preB.up, merged$FDR.q.val.iAMP.vs.preB.down, merged$FDR.q.val.PC.vs.preB.up, merged$FDR.q.val.PC.vs.preB.down, merged$FDR.q.val.ER.vs.preB.up, merged$FDR.q.val.ER.vs.preB.down, 
								merged$FDR.q.val.iAMP.vs.mature.up, merged$FDR.q.val.iAMP.vs.mature.down, merged$FDR.q.val.PC.vs.mature.up, merged$FDR.q.val.PC.vs.mature.down, merged$FDR.q.val.ER.vs.mature.up, merged$FDR.q.val.ER.vs.mature.down, 
								na.rm=T)

# NES = normalized enrichment score = value in heatmap
# stretch NES range for better color separation in heatmap

#merged$iAMP.vs.ER <- ifelse(!is.na(merged$FDR.q.val.iAMP.vs.ER.up) & merged$FDR.q.val.iAMP.vs.ER.up < 0.1, merged$NES.iAMP.vs.ER.up, ifelse(!is.na(merged$FDR.q.val.iAMP.vs.ER.down) & merged$FDR.q.val.iAMP.vs.ER.down < 0.1, merged$NES.iAMP.vs.ER.down, 0))
merged$iAMP.vs.ER <- ifelse(!is.na(merged$NES.iAMP.vs.ER.up), merged$NES.iAMP.vs.ER.up, merged$NES.iAMP.vs.ER.down)
merged$FDR.q.val.iAMP.vs.ER <- pmin(merged$FDR.q.val.iAMP.vs.ER.up, merged$FDR.q.val.iAMP.vs.ER.down, na.rm=T)
merged$iAMP.vs.ER[is.na(merged$iAMP.vs.ER)] <- 0
merged$iAMP.vs.ER <- ifelse(merged$iAMP.vs.ER > 0, merged$iAMP.vs.ER*merged$iAMP.vs.ER, -merged$iAMP.vs.ER*merged$iAMP.vs.ER)
#merged$iAMP.vs.PC <- ifelse(!is.na(merged$FDR.q.val.iAMP.vs.PC.up) & merged$FDR.q.val.iAMP.vs.PC.up < 0.1, merged$NES.iAMP.vs.PC.up, ifelse(!is.na(merged$FDR.q.val.iAMP.vs.PC.down) & merged$FDR.q.val.iAMP.vs.PC.down < 0.1, merged$NES.iAMP.vs.PC.down, 0))
merged$iAMP.vs.PC <- ifelse(!is.na(merged$NES.iAMP.vs.PC.up), merged$NES.iAMP.vs.PC.up, merged$NES.iAMP.vs.PC.down)
merged$FDR.q.val.iAMP.vs.PC <- pmin(merged$FDR.q.val.iAMP.vs.PC.up, merged$FDR.q.val.iAMP.vs.PC.down, na.rm=T)
merged$iAMP.vs.PC[is.na(merged$iAMP.vs.PC)] <- 0
merged$iAMP.vs.PC <- ifelse(merged$iAMP.vs.PC > 0, merged$iAMP.vs.PC*merged$iAMP.vs.PC, -merged$iAMP.vs.PC*merged$iAMP.vs.PC)
#merged$ER.vs.PC <- ifelse(!is.na(merged$FDR.q.val.ER.vs.PC.up) & merged$FDR.q.val.ER.vs.PC.up < 0.1, merged$NES.ER.vs.PC.up, ifelse(!is.na(merged$FDR.q.val.ER.vs.PC.down) & merged$FDR.q.val.ER.vs.PC.down < 0.1, merged$NES.ER.vs.PC.down, 0))
merged$ER.vs.PC <- ifelse(!is.na(merged$NES.ER.vs.PC.up), merged$NES.ER.vs.PC.up, merged$NES.ER.vs.PC.down)
merged$FDR.q.val.ER.vs.PC <- pmin(merged$FDR.q.val.ER.vs.PC.up, merged$FDR.q.val.ER.vs.PC.down, na.rm=T)
merged$ER.vs.PC[is.na(merged$ER.vs.PC)] <- 0
merged$ER.vs.PC <- ifelse(merged$ER.vs.PC > 0, merged$ER.vs.PC*merged$ER.vs.PC, -merged$ER.vs.PC*merged$ER.vs.PC)
merged$iAMP.vs.immature <- ifelse(!is.na(merged$NES.iAMP.vs.immature.up), merged$NES.iAMP.vs.immature.up, merged$NES.iAMP.vs.immature.down)
merged$FDR.q.val.iAMP.vs.immature <- pmin(merged$FDR.q.val.iAMP.vs.immature.up, merged$FDR.q.val.iAMP.vs.immature.down, na.rm=T)
merged$iAMP.vs.immature[is.na(merged$iAMP.vs.immature)] <- 0
merged$iAMP.vs.immature <- ifelse(merged$iAMP.vs.immature > 0, merged$iAMP.vs.immature*merged$iAMP.vs.immature, -merged$iAMP.vs.immature*merged$iAMP.vs.immature)
merged$PC.vs.immature <- ifelse(!is.na(merged$NES.PC.vs.immature.up), merged$NES.PC.vs.immature.up, merged$NES.PC.vs.immature.down)
merged$FDR.q.val.PC.vs.immature <- pmin(merged$FDR.q.val.PC.vs.immature.up, merged$FDR.q.val.PC.vs.immature.down, na.rm=T)
merged$PC.vs.immature[is.na(merged$PC.vs.immature)] <- 0
merged$PC.vs.immature <- ifelse(merged$PC.vs.immature > 0, merged$PC.vs.immature*merged$PC.vs.immature, -merged$PC.vs.immature*merged$PC.vs.immature)
merged$ER.vs.immature <- ifelse(!is.na(merged$NES.ER.vs.immature.up), merged$NES.ER.vs.immature.up, merged$NES.ER.vs.immature.down)
merged$FDR.q.val.ER.vs.immature <- pmin(merged$FDR.q.val.ER.vs.immature.up, merged$FDR.q.val.ER.vs.immature.down, na.rm=T)
merged$ER.vs.immature[is.na(merged$ER.vs.immature)] <- 0
merged$ER.vs.immature <- ifelse(merged$ER.vs.immature > 0, merged$ER.vs.immature*merged$ER.vs.immature, -merged$ER.vs.immature*merged$ER.vs.immature)
merged$iAMP.vs.preB <- ifelse(!is.na(merged$NES.iAMP.vs.preB.up), merged$NES.iAMP.vs.preB.up, merged$NES.iAMP.vs.preB.down)
merged$FDR.q.val.iAMP.vs.preB <- pmin(merged$FDR.q.val.iAMP.vs.preB.up, merged$FDR.q.val.iAMP.vs.preB.down, na.rm=T)
merged$iAMP.vs.preB[is.na(merged$iAMP.vs.preB)] <- 0
merged$iAMP.vs.preB <- ifelse(merged$iAMP.vs.preB > 0, merged$iAMP.vs.preB*merged$iAMP.vs.preB, -merged$iAMP.vs.preB*merged$iAMP.vs.preB)
merged$PC.vs.preB <- ifelse(!is.na(merged$NES.PC.vs.preB.up), merged$NES.PC.vs.preB.up, merged$NES.PC.vs.preB.down)
merged$FDR.q.val.PC.vs.preB <- pmin(merged$FDR.q.val.PC.vs.preB.up, merged$FDR.q.val.PC.vs.preB.down, na.rm=T)
merged$PC.vs.preB[is.na(merged$PC.vs.preB)] <- 0
merged$PC.vs.preB <- ifelse(merged$PC.vs.preB > 0, merged$PC.vs.preB*merged$PC.vs.preB, -merged$PC.vs.preB*merged$PC.vs.preB)
merged$ER.vs.preB <- ifelse(!is.na(merged$NES.ER.vs.preB.up), merged$NES.ER.vs.preB.up, merged$NES.ER.vs.preB.down)
merged$FDR.q.val.ER.vs.preB <- pmin(merged$FDR.q.val.ER.vs.preB.up, merged$FDR.q.val.ER.vs.preB.down, na.rm=T)
merged$ER.vs.preB[is.na(merged$ER.vs.preB)] <- 0
merged$ER.vs.preB <- ifelse(merged$ER.vs.preB > 0, merged$ER.vs.preB*merged$ER.vs.preB, -merged$ER.vs.preB*merged$ER.vs.preB)
merged$iAMP.vs.mature <- ifelse(!is.na(merged$NES.iAMP.vs.mature.up), merged$NES.iAMP.vs.mature.up, merged$NES.iAMP.vs.mature.down)
merged$FDR.q.val.iAMP.vs.mature <- pmin(merged$FDR.q.val.iAMP.vs.mature.up, merged$FDR.q.val.iAMP.vs.mature.down, na.rm=T)
merged$iAMP.vs.mature[is.na(merged$iAMP.vs.mature)] <- 0
merged$iAMP.vs.mature <- ifelse(merged$iAMP.vs.mature > 0, merged$iAMP.vs.mature*merged$iAMP.vs.mature, -merged$iAMP.vs.mature*merged$iAMP.vs.mature)
merged$PC.vs.mature <- ifelse(!is.na(merged$NES.PC.vs.mature.up), merged$NES.PC.vs.mature.up, merged$NES.PC.vs.mature.down)
merged$FDR.q.val.PC.vs.mature <- pmin(merged$FDR.q.val.PC.vs.mature.up, merged$FDR.q.val.PC.vs.mature.down, na.rm=T)
merged$PC.vs.mature[is.na(merged$PC.vs.mature)] <- 0
merged$PC.vs.mature <- ifelse(merged$PC.vs.mature > 0, merged$PC.vs.mature*merged$PC.vs.mature, -merged$PC.vs.mature*merged$PC.vs.mature)
merged$ER.vs.mature <- ifelse(!is.na(merged$NES.ER.vs.mature.up), merged$NES.ER.vs.mature.up, merged$NES.ER.vs.mature.down)
merged$FDR.q.val.ER.vs.mature <- pmin(merged$FDR.q.val.ER.vs.mature.up, merged$FDR.q.val.ER.vs.mature.down, na.rm=T)
merged$ER.vs.mature[is.na(merged$ER.vs.mature)] <- 0
merged$ER.vs.mature <- ifelse(merged$ER.vs.mature > 0, merged$ER.vs.mature*merged$ER.vs.mature, -merged$ER.vs.mature*merged$ER.vs.mature)

cols.values <- c("iAMP.vs.ER", "iAMP.vs.PC", "ER.vs.PC", "iAMP.vs.immature", "PC.vs.immature", "ER.vs.immature", "iAMP.vs.preB", "PC.vs.preB", "ER.vs.preB", "iAMP.vs.mature", "PC.vs.mature", "ER.vs.mature")

plot.heatmap <- function(data, cexCol=0.9, cexRow=0.9, sigLevel=0.02, minNES=NA, title="") {
	if (!is.na(minNES)) {
		nes.best <- apply(data[,c("iAMP.vs.ER", "iAMP.vs.PC", "ER.vs.PC", "iAMP.vs.immature", "PC.vs.immature", "ER.vs.immature", "iAMP.vs.preB", "PC.vs.preB", "ER.vs.preB", "iAMP.vs.mature", "PC.vs.mature", "ER.vs.mature")], 1, function (x) { max(abs(x))} )
		data <- data[nes.best >= minNES,]
	}
	data <- data[data$FDR.best<sigLevel,cols.values]
	hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
	fdrs <- as.matrix(merged[rownames(data),c("FDR.q.val.iAMP.vs.ER", "FDR.q.val.iAMP.vs.PC", "FDR.q.val.ER.vs.PC", "FDR.q.val.iAMP.vs.immature", "FDR.q.val.PC.vs.immature", "FDR.q.val.ER.vs.immature", "FDR.q.val.iAMP.vs.preB", "FDR.q.val.PC.vs.preB", "FDR.q.val.ER.vs.preB", "FDR.q.val.iAMP.vs.mature", "FDR.q.val.PC.vs.mature", "FDR.q.val.ER.vs.mature")])
	fdrs[is.na(fdrs)|fdrs>sigLevel] <- NA
	fdrs[!is.na(fdrs)&fdrs<=sigLevel] <- "*"
	heatmap.2(as.matrix(data), Colv=F, Rowv=T, dendrogram="none", trace="none", col=rev(hmcol), margin=c(10, 25), cexCol=cexCol, cexRow=cexRow, keysize=0.7, 
			  colsep=seq(1:ncol(data)), rowsep=seq(1:nrow(data)), sepcolor="grey92", sepwidth=c(0.005,0.005),
			  cellnote=fdrs, notecol='white',
			  main=sprintf("%s (FDR <= %.4g)", title, sigLevel))
	
}

# PATHWAYS
pdf("/mnt/projects/iamp/results/gsea-heatmap.pathways.pdf", height=15, width=10)
geneset <- grep("^(REACTOME|KEGG|PID_)", rownames(merged), perl=T, value=TRUE)
genesets <- geneset
plot.heatmap(merged[geneset,], cexRow=0.7, sigLevel=0.001, title="Pathways")
dev.off()

# UP
pdf("/mnt/projects/iamp/results/gsea-heatmap.up.pdf", height=15, width=10)
geneset <- grep("_UP$", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, sigLevel=0.001, minNES=7, title="UP genesets")
dev.off()

# DOWN
pdf("/mnt/projects/iamp/results/gsea-heatmap.dn.pdf", height=15, width=10)
geneset <- grep("_DN$", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, sigLevel=0.001, minNES=7, title="DOWN genesets")
dev.off()

# MIR TARGETS
pdf("/mnt/projects/iamp/results/gsea-heatmap.mir.pdf", height=15, width=10)
geneset <- grep("MIR-", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], sigLevel=0.05, title="MiR targets")
dev.off()

# CHROMOSOME POSITION
pdf("/mnt/projects/iamp/results/gsea-heatmap.chr.pdf", height=15, width=10)
geneset <- grep("^CHR[\\dXY]", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], title="Chromosome location")
dev.off()

# SIGNALING
pdf("/mnt/projects/iamp/results/gsea-heatmap.signaling.pdf", height=15, width=10)
geneset <- grep("SIGNALING", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], title="Signaling")
dev.off()

# PROMOTER MOTIF
pdf("/mnt/projects/iamp/results/gsea-heatmap.promoter-motif.pdf", height=15, width=10)
geneset <- grep("V\\$", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, sigLevel=0.01, title="Motif in promoter region")
dev.off()

# MORF expression compendium
pdf("/mnt/projects/iamp/results/gsea-heatmap.morf.pdf", height=15, width=10)
geneset <- grep("^MORF_", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, title="MORF neighborhood")
dev.off()

# MODULES
pdf("/mnt/projects/iamp/results/gsea-heatmap.module.pdf", height=15, width=10)
geneset <- grep("^MODULE_", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, sigLevel=0.001, title="Modules")
dev.off()

# OTHER
pdf("/mnt/projects/iamp/results/gsea-heatmap.other.pdf", height=15, width=10)
geneset <- rownames(merged)[!rownames(merged) %in% genesets]
plot.heatmap(merged[geneset,], cexRow=0.7, sigLevel=0.0001, minNES=7, title="Other")
dev.off()

