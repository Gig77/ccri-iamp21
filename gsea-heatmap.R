library(RColorBrewer)
library("gplots")

iAMP.vs.ER.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsER-gsea/enrichedUp.csv")
iAMP.vs.ER.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsER-gsea/enrichedDown.csv")
iAMP.vs.DS.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsDS-gsea/enrichedUp.csv")
iAMP.vs.DS.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsDS-gsea/enrichedDown.csv")
ER.vs.DS.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsDS-gsea/enrichedUp.csv")
ER.vs.DS.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsDS-gsea/enrichedDown.csv")
iAMP.vs.noniAMP.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsNoniAMP-gsea/enrichedUp.csv")
iAMP.vs.noniAMP.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsNoniAMP-gsea/enrichedDown.csv")
DS.vs.nonDS.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsNonDS-gsea/enrichedUp.csv")
DS.vs.nonDS.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsNonDS-gsea/enrichedDown.csv")
ER.vs.nonER.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsNonER-gsea/enrichedUp.csv")
ER.vs.nonER.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsNonER-gsea/enrichedDown.csv")
iAMP.vs.immature.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsImmature-gsea/enrichedUp.csv")
iAMP.vs.immature.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsImmature-gsea/enrichedDown.csv")
DS.vs.immature.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsImmature-gsea/enrichedUp.csv")
DS.vs.immature.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsImmature-gsea/enrichedDown.csv")
ER.vs.immature.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsImmature-gsea/enrichedUp.csv")
ER.vs.immature.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsImmature-gsea/enrichedDown.csv")
iAMP.vs.preB.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsPreB-gsea/enrichedUp.csv")
iAMP.vs.preB.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsPreB-gsea/enrichedDown.csv")
DS.vs.preB.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsPreB-gsea/enrichedUp.csv")
DS.vs.preB.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsPreB-gsea/enrichedDown.csv")
ER.vs.preB.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsPreB-gsea/enrichedUp.csv")
ER.vs.preB.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsPreB-gsea/enrichedDown.csv")
iAMP.vs.mature.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsMature-gsea/enrichedUp.csv")
iAMP.vs.mature.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_iAMPvsMature-gsea/enrichedDown.csv")
DS.vs.mature.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsMature-gsea/enrichedUp.csv")
DS.vs.mature.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_DSvsMature-gsea/enrichedDown.csv")
ER.vs.mature.up <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsMature-gsea/enrichedUp.csv")
ER.vs.mature.down <- read.delim("/mnt/projects/iamp/results/anduril/execute/gsea_ERvsMature-gsea/enrichedDown.csv")

cols.keep <- c(1, 5:11)
cols.merge <- c(1)
merged <- merge(iAMP.vs.ER.up[,cols.keep], iAMP.vs.ER.down[,cols.keep], by=cols.merge, suffixes=c(".iAMP.vs.ER.up", ".iAMP.vs.ER.down"), all=T)
merged <- merge(merged, iAMP.vs.DS.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.DS.up")
merged <- merge(merged, iAMP.vs.DS.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.DS.down")
merged <- merge(merged, ER.vs.DS.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.DS.up")
merged <- merge(merged, ER.vs.DS.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.DS.down")
merged <- merge(merged, iAMP.vs.noniAMP.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.noniAMP.up")
merged <- merge(merged, iAMP.vs.noniAMP.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.noniAMP.down")
merged <- merge(merged, DS.vs.nonDS.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.nonDS.up")
merged <- merge(merged, DS.vs.nonDS.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.nonDS.down")
merged <- merge(merged, ER.vs.nonER.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.nonER.up")
merged <- merge(merged, ER.vs.nonER.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.nonER.down")
merged <- merge(merged, iAMP.vs.immature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.immature.up")
merged <- merge(merged, iAMP.vs.immature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.immature.down")
merged <- merge(merged, DS.vs.immature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.immature.up")
merged <- merge(merged, DS.vs.immature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.immature.down")
merged <- merge(merged, ER.vs.immature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.immature.up")
merged <- merge(merged, ER.vs.immature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.immature.down")
merged <- merge(merged, iAMP.vs.preB.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.preB.up")
merged <- merge(merged, iAMP.vs.preB.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.preB.down")
merged <- merge(merged, DS.vs.preB.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.preB.up")
merged <- merge(merged, DS.vs.preB.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.preB.down")
merged <- merge(merged, ER.vs.preB.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.preB.up")
merged <- merge(merged, ER.vs.preB.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.preB.down")
merged <- merge(merged, iAMP.vs.mature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.mature.up")
merged <- merge(merged, iAMP.vs.mature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".iAMP.vs.mature.down")
merged <- merge(merged, DS.vs.mature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.mature.up")
merged <- merge(merged, DS.vs.mature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".DS.vs.mature.down")
merged <- merge(merged, ER.vs.mature.up[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.mature.up")
merged <- merge(merged, ER.vs.mature.down[,cols.keep], by=cols.merge, all=T) ; names(merged)[(ncol(merged)-6):ncol(merged)] <- paste0(names(merged)[(ncol(merged)-6):ncol(merged)], ".ER.vs.mature.down")
rownames(merged) <- merged$NAME

# get top FDR for gene set
#merged$FDR.best <- mapply(pmin, merged$FDR.q.val.iAMP.vs.ER.up, merged$FDR.q.val.iAMP.vs.ER.down, merged$FDR.q.val.iAMP.vs.DS.up, merged$FDR.q.val.iAMP.vs.DS.down, merged$FDR.q.val.ER.vs.DS.up, merged$FDR.q.val.ER.vs.DS.down, 
#		                        merged$FDR.q.val.iAMP.vs.immature.up, merged$FDR.q.val.iAMP.vs.immature.down, merged$FDR.q.val.DS.vs.immature.up, merged$FDR.q.val.DS.vs.immature.down, merged$FDR.q.val.ER.vs.immature.up, merged$FDR.q.val.ER.vs.immature.down, 
#								merged$FDR.q.val.iAMP.vs.preB.up, merged$FDR.q.val.iAMP.vs.preB.down, merged$FDR.q.val.DS.vs.preB.up, merged$FDR.q.val.DS.vs.preB.down, merged$FDR.q.val.ER.vs.preB.up, merged$FDR.q.val.ER.vs.preB.down, 
#								merged$FDR.q.val.iAMP.vs.mature.up, merged$FDR.q.val.iAMP.vs.mature.down, merged$FDR.q.val.DS.vs.mature.up, merged$FDR.q.val.DS.vs.mature.down, merged$FDR.q.val.ER.vs.mature.up, merged$FDR.q.val.ER.vs.mature.down, 
#								na.rm=T)

merged$FDR.best <- mapply(pmin, merged$FDR.q.val.iAMP.vs.ER.up,      merged$FDR.q.val.iAMP.vs.ER.down, 
                                merged$FDR.q.val.iAMP.vs.DS.up,      merged$FDR.q.val.iAMP.vs.DS.down, 
                                merged$FDR.q.val.ER.vs.DS.up,        merged$FDR.q.val.ER.vs.DS.down, 
                                merged$FDR.q.val.iAMP.vs.noniAMP.up, merged$FDR.q.val.iAMP.vs.noniAMP.down, 
                                merged$FDR.q.val.DS.vs.nonDS.up, merged$FDR.q.val.DS.vs.nonDS.down, 
                                merged$FDR.q.val.ER.vs.nonER.up, merged$FDR.q.val.ER.vs.nonER.down, 
                          na.rm=T)
merged <- merged[!is.na(merged$FDR.best),]

# NES = normalized enrichment score = value in heatmap
# stretch NES range for better color separation in heatmap

#merged$iAMP.vs.ER <- ifelse(!is.na(merged$FDR.q.val.iAMP.vs.ER.up) & merged$FDR.q.val.iAMP.vs.ER.up < 0.1, merged$NES.iAMP.vs.ER.up, ifelse(!is.na(merged$FDR.q.val.iAMP.vs.ER.down) & merged$FDR.q.val.iAMP.vs.ER.down < 0.1, merged$NES.iAMP.vs.ER.down, 0))
merged$iAMP.vs.ER <- ifelse(!is.na(merged$NES.iAMP.vs.ER.up), merged$NES.iAMP.vs.ER.up, merged$NES.iAMP.vs.ER.down)
merged$FDR.q.val.iAMP.vs.ER <- pmin(merged$FDR.q.val.iAMP.vs.ER.up, merged$FDR.q.val.iAMP.vs.ER.down, na.rm=T)
merged$iAMP.vs.ER[is.na(merged$iAMP.vs.ER)] <- 0
merged$iAMP.vs.ER <- ifelse(merged$iAMP.vs.ER > 0, merged$iAMP.vs.ER*merged$iAMP.vs.ER, -merged$iAMP.vs.ER*merged$iAMP.vs.ER)

#merged$iAMP.vs.DS <- ifelse(!is.na(merged$FDR.q.val.iAMP.vs.DS.up) & merged$FDR.q.val.iAMP.vs.DS.up < 0.1, merged$NES.iAMP.vs.DS.up, ifelse(!is.na(merged$FDR.q.val.iAMP.vs.DS.down) & merged$FDR.q.val.iAMP.vs.DS.down < 0.1, merged$NES.iAMP.vs.DS.down, 0))
merged$iAMP.vs.DS <- ifelse(!is.na(merged$NES.iAMP.vs.DS.up), merged$NES.iAMP.vs.DS.up, merged$NES.iAMP.vs.DS.down)
merged$FDR.q.val.iAMP.vs.DS <- pmin(merged$FDR.q.val.iAMP.vs.DS.up, merged$FDR.q.val.iAMP.vs.DS.down, na.rm=T)
merged$iAMP.vs.DS[is.na(merged$iAMP.vs.DS)] <- 0
merged$iAMP.vs.DS <- ifelse(merged$iAMP.vs.DS > 0, merged$iAMP.vs.DS*merged$iAMP.vs.DS, -merged$iAMP.vs.DS*merged$iAMP.vs.DS)

#merged$ER.vs.DS <- ifelse(!is.na(merged$FDR.q.val.ER.vs.DS.up) & merged$FDR.q.val.ER.vs.DS.up < 0.1, merged$NES.ER.vs.DS.up, ifelse(!is.na(merged$FDR.q.val.ER.vs.DS.down) & merged$FDR.q.val.ER.vs.DS.down < 0.1, merged$NES.ER.vs.DS.down, 0))
merged$ER.vs.DS <- ifelse(!is.na(merged$NES.ER.vs.DS.up), merged$NES.ER.vs.DS.up, merged$NES.ER.vs.DS.down)
merged$FDR.q.val.ER.vs.DS <- pmin(merged$FDR.q.val.ER.vs.DS.up, merged$FDR.q.val.ER.vs.DS.down, na.rm=T)
merged$ER.vs.DS[is.na(merged$ER.vs.DS)] <- 0
merged$ER.vs.DS <- ifelse(merged$ER.vs.DS > 0, merged$ER.vs.DS*merged$ER.vs.DS, -merged$ER.vs.DS*merged$ER.vs.DS)

merged$iAMP.vs.noniAMP <- ifelse(!is.na(merged$NES.iAMP.vs.noniAMP.up), merged$NES.iAMP.vs.noniAMP.up, merged$NES.iAMP.vs.noniAMP.down)
merged$FDR.q.val.iAMP.vs.noniAMP <- pmin(merged$FDR.q.val.iAMP.vs.noniAMP.up, merged$FDR.q.val.iAMP.vs.noniAMP.down, na.rm=T)
merged$iAMP.vs.noniAMP[is.na(merged$iAMP.vs.noniAMP)] <- 0
merged$iAMP.vs.noniAMP <- ifelse(merged$iAMP.vs.noniAMP > 0, merged$iAMP.vs.noniAMP*merged$iAMP.vs.noniAMP, -merged$iAMP.vs.noniAMP*merged$iAMP.vs.noniAMP)

merged$DS.vs.nonDS <- ifelse(!is.na(merged$NES.DS.vs.nonDS.up), merged$NES.DS.vs.nonDS.up, merged$NES.DS.vs.nonDS.down)
merged$FDR.q.val.DS.vs.nonDS <- pmin(merged$FDR.q.val.DS.vs.nonDS.up, merged$FDR.q.val.DS.vs.nonDS.down, na.rm=T)
merged$DS.vs.nonDS[is.na(merged$DS.vs.nonDS)] <- 0
merged$DS.vs.nonDS <- ifelse(merged$DS.vs.nonDS > 0, merged$DS.vs.nonDS*merged$DS.vs.nonDS, -merged$DS.vs.nonDS*merged$DS.vs.nonDS)

merged$ER.vs.nonER <- ifelse(!is.na(merged$NES.ER.vs.nonER.up), merged$NES.ER.vs.nonER.up, merged$NES.ER.vs.nonER.down)
merged$FDR.q.val.ER.vs.nonER <- pmin(merged$FDR.q.val.ER.vs.nonER.up, merged$FDR.q.val.ER.vs.nonER.down, na.rm=T)
merged$ER.vs.nonER[is.na(merged$ER.vs.nonER)] <- 0
merged$ER.vs.nonER <- ifelse(merged$ER.vs.nonER > 0, merged$ER.vs.nonER*merged$ER.vs.nonER, -merged$ER.vs.nonER*merged$ER.vs.nonER)

merged$iAMP.vs.immature <- ifelse(!is.na(merged$NES.iAMP.vs.immature.up), merged$NES.iAMP.vs.immature.up, merged$NES.iAMP.vs.immature.down)
merged$FDR.q.val.iAMP.vs.immature <- pmin(merged$FDR.q.val.iAMP.vs.immature.up, merged$FDR.q.val.iAMP.vs.immature.down, na.rm=T)
merged$iAMP.vs.immature[is.na(merged$iAMP.vs.immature)] <- 0
merged$iAMP.vs.immature <- ifelse(merged$iAMP.vs.immature > 0, merged$iAMP.vs.immature*merged$iAMP.vs.immature, -merged$iAMP.vs.immature*merged$iAMP.vs.immature)

merged$DS.vs.immature <- ifelse(!is.na(merged$NES.DS.vs.immature.up), merged$NES.DS.vs.immature.up, merged$NES.DS.vs.immature.down)
merged$FDR.q.val.DS.vs.immature <- pmin(merged$FDR.q.val.DS.vs.immature.up, merged$FDR.q.val.DS.vs.immature.down, na.rm=T)
merged$DS.vs.immature[is.na(merged$DS.vs.immature)] <- 0
merged$DS.vs.immature <- ifelse(merged$DS.vs.immature > 0, merged$DS.vs.immature*merged$DS.vs.immature, -merged$DS.vs.immature*merged$DS.vs.immature)

merged$ER.vs.immature <- ifelse(!is.na(merged$NES.ER.vs.immature.up), merged$NES.ER.vs.immature.up, merged$NES.ER.vs.immature.down)
merged$FDR.q.val.ER.vs.immature <- pmin(merged$FDR.q.val.ER.vs.immature.up, merged$FDR.q.val.ER.vs.immature.down, na.rm=T)
merged$ER.vs.immature[is.na(merged$ER.vs.immature)] <- 0
merged$ER.vs.immature <- ifelse(merged$ER.vs.immature > 0, merged$ER.vs.immature*merged$ER.vs.immature, -merged$ER.vs.immature*merged$ER.vs.immature)

merged$iAMP.vs.preB <- ifelse(!is.na(merged$NES.iAMP.vs.preB.up), merged$NES.iAMP.vs.preB.up, merged$NES.iAMP.vs.preB.down)
merged$FDR.q.val.iAMP.vs.preB <- pmin(merged$FDR.q.val.iAMP.vs.preB.up, merged$FDR.q.val.iAMP.vs.preB.down, na.rm=T)
merged$iAMP.vs.preB[is.na(merged$iAMP.vs.preB)] <- 0
merged$iAMP.vs.preB <- ifelse(merged$iAMP.vs.preB > 0, merged$iAMP.vs.preB*merged$iAMP.vs.preB, -merged$iAMP.vs.preB*merged$iAMP.vs.preB)

merged$DS.vs.preB <- ifelse(!is.na(merged$NES.DS.vs.preB.up), merged$NES.DS.vs.preB.up, merged$NES.DS.vs.preB.down)
merged$FDR.q.val.DS.vs.preB <- pmin(merged$FDR.q.val.DS.vs.preB.up, merged$FDR.q.val.DS.vs.preB.down, na.rm=T)
merged$DS.vs.preB[is.na(merged$DS.vs.preB)] <- 0
merged$DS.vs.preB <- ifelse(merged$DS.vs.preB > 0, merged$DS.vs.preB*merged$DS.vs.preB, -merged$DS.vs.preB*merged$DS.vs.preB)

merged$ER.vs.preB <- ifelse(!is.na(merged$NES.ER.vs.preB.up), merged$NES.ER.vs.preB.up, merged$NES.ER.vs.preB.down)
merged$FDR.q.val.ER.vs.preB <- pmin(merged$FDR.q.val.ER.vs.preB.up, merged$FDR.q.val.ER.vs.preB.down, na.rm=T)
merged$ER.vs.preB[is.na(merged$ER.vs.preB)] <- 0
merged$ER.vs.preB <- ifelse(merged$ER.vs.preB > 0, merged$ER.vs.preB*merged$ER.vs.preB, -merged$ER.vs.preB*merged$ER.vs.preB)

merged$iAMP.vs.mature <- ifelse(!is.na(merged$NES.iAMP.vs.mature.up), merged$NES.iAMP.vs.mature.up, merged$NES.iAMP.vs.mature.down)
merged$FDR.q.val.iAMP.vs.mature <- pmin(merged$FDR.q.val.iAMP.vs.mature.up, merged$FDR.q.val.iAMP.vs.mature.down, na.rm=T)
merged$iAMP.vs.mature[is.na(merged$iAMP.vs.mature)] <- 0
merged$iAMP.vs.mature <- ifelse(merged$iAMP.vs.mature > 0, merged$iAMP.vs.mature*merged$iAMP.vs.mature, -merged$iAMP.vs.mature*merged$iAMP.vs.mature)

merged$DS.vs.mature <- ifelse(!is.na(merged$NES.DS.vs.mature.up), merged$NES.DS.vs.mature.up, merged$NES.DS.vs.mature.down)
merged$FDR.q.val.DS.vs.mature <- pmin(merged$FDR.q.val.DS.vs.mature.up, merged$FDR.q.val.DS.vs.mature.down, na.rm=T)
merged$DS.vs.mature[is.na(merged$DS.vs.mature)] <- 0
merged$DS.vs.mature <- ifelse(merged$DS.vs.mature > 0, merged$DS.vs.mature*merged$DS.vs.mature, -merged$DS.vs.mature*merged$DS.vs.mature)

merged$ER.vs.mature <- ifelse(!is.na(merged$NES.ER.vs.mature.up), merged$NES.ER.vs.mature.up, merged$NES.ER.vs.mature.down)
merged$FDR.q.val.ER.vs.mature <- pmin(merged$FDR.q.val.ER.vs.mature.up, merged$FDR.q.val.ER.vs.mature.down, na.rm=T)
merged$ER.vs.mature[is.na(merged$ER.vs.mature)] <- 0
merged$ER.vs.mature <- ifelse(merged$ER.vs.mature > 0, merged$ER.vs.mature*merged$ER.vs.mature, -merged$ER.vs.mature*merged$ER.vs.mature)

cols.values <- c("iAMP.vs.noniAMP", "DS.vs.nonDS", "ER.vs.nonER", "iAMP.vs.ER", "iAMP.vs.DS", "ER.vs.DS", "iAMP.vs.immature", "DS.vs.immature", "ER.vs.immature", "iAMP.vs.preB", "DS.vs.preB", "ER.vs.preB", "iAMP.vs.mature", "DS.vs.mature", "ER.vs.mature")

plot.heatmap <- function(data, cexCol=0.9, cexRow=0.9, sigLevel=0.02, hsigLevel=1e-5, minNES=NA, title="") {
	if (!is.na(minNES)) {
		nes.best <- apply(data[,c("iAMP.vs.noniAMP", "DS.vs.nonDS", "ER.vs.nonER", "iAMP.vs.ER", "iAMP.vs.DS", "ER.vs.DS")], 1, function (x) { max(abs(x))} )
		data <- data[nes.best >= minNES,]
	}
	data <- data[data$FDR.best<sigLevel,cols.values]
	hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
	fdrs <- as.matrix(merged[rownames(data),c("FDR.q.val.iAMP.vs.noniAMP", "FDR.q.val.DS.vs.nonDS", "FDR.q.val.ER.vs.nonER", "FDR.q.val.iAMP.vs.ER", "FDR.q.val.iAMP.vs.DS", "FDR.q.val.ER.vs.DS", "FDR.q.val.iAMP.vs.immature", "FDR.q.val.DS.vs.immature", "FDR.q.val.ER.vs.immature", "FDR.q.val.iAMP.vs.preB", "FDR.q.val.DS.vs.preB", "FDR.q.val.ER.vs.preB", "FDR.q.val.iAMP.vs.mature", "FDR.q.val.DS.vs.mature", "FDR.q.val.ER.vs.mature")])
	sig <- fdrs ; sig[,] <- NA
	sig[fdrs <= sigLevel & fdrs > hsigLevel] <- "*"
	sig[fdrs <= hsigLevel] <- "**"
	heatmap.2(as.matrix(data), Colv=F, Rowv=T, dendrogram="row", trace="none", col=rev(hmcol), margin=c(10, 25), cexCol=cexCol, cexRow=cexRow, keysize=0.7, 
			  colsep=seq(1:ncol(data)), rowsep=seq(1:nrow(data)), sepcolor="grey92", sepwidth=c(0.005,0.005),
			  cellnote=sig, notecol='white',
			  main=sprintf("%s\n(* FDR <= %.4g, ** FDR <= %.4g%s)", title, sigLevel, hsigLevel, ifelse(!is.na(minNES), sprintf(",\nabs(NES) >= %.1f", minNES), "")),
			  key.title="NES")
	
}

# PATHWAYS
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_pathways.pdf", height=15, width=10)
geneset <- grep("^(REACTOME|KEGG|PID_|BIOCARTA)", rownames(merged), perl=T, value=TRUE)
genesets <- geneset
plot.heatmap(merged[geneset,], cexRow=0.7, sigLevel=0.001, hsigLevel=1e-10, title="MSigDB 5.0 Pathways")
dev.off()

# CHEMICAL PERTURBATION
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_chemical-perturbation.pdf", height=15, width=10)
geneset <- grep("(_UP|_DN)$", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.7, sigLevel=1e-3, hsigLevel=1e-10, minNES=5.5, title="MSigDB 5.0 Chemical Perturbation")
dev.off()

# MIR TARGETS
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_mir.pdf", height=15, width=10)
geneset <- grep("MIR-", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], sigLevel=0.1, hsigLevel=1e-2, title="MSigDB 5.0 MiR Targets")
dev.off()

# CUSTOM MIR TARGET GENE SETS
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_custom.pdf", height=8, width=10)
geneset <- grep("(PREDICTED|VALIDATED|HENNING)", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], sigLevel=0.3, hsigLevel=0.05, title="Custom MiR Target Gene Sets")
dev.off()

# CHROMOSOME POSITION
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_chr.pdf", height=15, width=10)
geneset <- grep("^CHR[\\dXY]", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], title="Chromosome Location", cexRow=1, sigLevel=1e-3, hsigLevel=1e-8)
dev.off()

# SIGNALING
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_signaling.pdf", height=15, width=10)
geneset <- grep("SIGNALING", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], title="MSigDB 5.0 Signaling", sigLevel=0.01, hsigLevel=1e-10)
dev.off()

# TF TARGETS
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_TF-targets.pdf", height=15, width=10)
geneset <- grep("V\\$", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.8, sigLevel=0.01, title="MSigDB 5.0 TF target gene sets")
dev.off()

# Cancer gene neighborhoods
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_neighborhood.pdf", height=15, width=10)
geneset <- grep("^(CAR_|GCM_|GNF2_|MORF_)", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, sigLevel=1e-10, hsigLevel=1e-10, title="MSigDB 5.0 Cancer Gene Neighborhoods")
dev.off()

# MODULES
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_module.pdf", height=15, width=10)
geneset <- grep("^MODULE_", rownames(merged), perl=T, value=TRUE)
genesets <- c(genesets, geneset)
plot.heatmap(merged[geneset,], cexRow=0.5, sigLevel=0.001, title="MSigDB 5.0 Modules")
dev.off()

# OTHER
pdf("/mnt/projects/iamp/results/gsea/msigdb5/gsea-heatmap_other.pdf", height=15, width=10)
geneset <- rownames(merged)[!rownames(merged) %in% genesets]
plot.heatmap(merged[geneset,], cexRow=0.7, sigLevel=1e-5, hsigLevel=1e-10, title="MSigDB 5.0 Other Gene Sets")
dev.off()

