validated <- read.delim("/mnt/projects/iamp/data/miRecords/miRecords.validatedTargets.txt", stringsAsFactors = FALSE)
predicted <- read.delim("/mnt/projects/iamp/data/miRecords/miRecords.predictedTargets.txt", stringsAsFactors = FALSE)

miR125.val <- sort(unique(gdata::trim(validated$symbol[grepl("miR-125[ab]", validated$miRNA_mature_ID)]))) # 64 targets
miRlet7c.val <- sort(unique(gdata::trim(validated$symbol[grepl("let-7c", validated$miRNA_mature_ID)]))) # 7 targets
miRlet7a.val <- sort(unique(gdata::trim(validated$symbol[grepl("let-7a", validated$miRNA_mature_ID)]))) # 23 targets
miRlet7e.val <- sort(unique(gdata::trim(validated$symbol[grepl("let-7e", validated$miRNA_mature_ID)]))) # 5 targets
miR99.val <- sort(unique(gdata::trim(validated$symbol[grepl("miR-99[ab]", validated$miRNA_mature_ID)]))) # 4 targets
miR155.val <- sort(unique(gdata::trim(validated$symbol[grepl("miR-155", validated$miRNA_mature_ID)]))) # 24 targets
miR100.val <- sort(unique(gdata::trim(validated$symbol[grepl("miR-100", validated$miRNA_mature_ID)]))) # 4 targets
miR802.val <- sort(unique(gdata::trim(validated$symbol[grepl("miR-802", validated$miRNA_mature_ID)]))) # 1 targets

miR100.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 3 & predicted$miRNA.ID %in% c("hsa-miR-100")]))) # 102 targets
miR125a.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 5 & predicted$miRNA.ID %in% c("hsa-miR-125a-3p", "hsa-miR-125a-5p")]))) # 141 targets
miR125b.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 5 & predicted$miRNA.ID %in% c("hsa-miR-125b")]))) # 134 targets
miR125.pred <- sort(unique(c(miR125a.pred, miR125b.pred))) # 196 targets
miR155.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 5 & predicted$miRNA.ID %in% c("hsa-miR-155")]))) # 137 targets
miR802.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 5 & predicted$miRNA.ID %in% c("hsa-miR-802")]))) # 218 targets
miR99a.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 3 & predicted$miRNA.ID %in% c("hsa-miR-99a")]))) # 271 targets
miR99b.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 3 & predicted$miRNA.ID %in% c("hsa-miR-99b")]))) # 204 targets
miR99.pred <- sort(unique(c(miR99a.pred, miR99b.pred))) # 380 targets
miRlet7c.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 5 & predicted$miRNA.ID %in% c("hsa-let-7c")]))) # 389 targets
miRlet7a.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 5 & predicted$miRNA.ID %in% c("hsa-let-7a")]))) # 390 targets
miRlet7e.pred <- sort(unique(gdata::trim(predicted$Symbol[predicted$Symbol != "" & predicted$databases >= 3 & predicted$miRNA.ID %in% c("hsa-let-7e")]))) # 1940 targets

# hennings excel sheet
predicted.henning <- read.delim("/mnt/projects/iamp/data/henning/miRRECORDS_let7-99-125_targets_2012_Abfrage.txt", stringsAsFactors = FALSE, na.strings = c("", "NA"))
miR125.henning <- sort(unique(predicted.henning$Symbol[!is.na(predicted.henning$X125.databases)])) # 2225 targets
miRlet7c.henning <- sort(unique(predicted.henning$Symbol[!is.na(predicted.henning$let7c.databases)])) # 720 targets
miR99.henning <- sort(unique(predicted.henning$Symbol[!is.na(predicted.henning$X99.databases)])) # 326 targets

nr <- max(length(miR125.val), length(miRlet7c.val), length(miRlet7a.val), length(miRlet7e.val), length(miR99.val), length(miR155.val), length(miR100.val), 
          length(miR802.val), length(miR125a.pred), length(miR125b.pred), length(miR125.pred), length(miRlet7c.pred), length(miRlet7a.pred), length(miRlet7e.pred),
          length(miR99a.pred), length(miR99b.pred), length(miR99.pred), length(miR155.pred), length(miR100.pred), length(miR802.pred),
          length(miR125.henning), length(miRlet7c.henning ), length(miR99.henning))

gmx <- data.frame(
  'mir_125_validated' = c(miR125.val, rep(NA, nr-length(miR125.val))),
  'let_7c_validated' = c(miRlet7c.val, rep(NA, nr-length(miRlet7c.val))),
  'let_7a_validated' = c(miRlet7a.val, rep(NA, nr-length(miRlet7a.val))),
  'let_7e_validated' = c(miRlet7e.val, rep(NA, nr-length(miRlet7e.val))),
  'mir_99_validated' = c(miR99.val, rep(NA, nr-length(miR99.val))),
  'mir_155_validated' = c(miR155.val, rep(NA, nr-length(miR155.val))),
  'mir_100_validated' = c(miR100.val, rep(NA, nr-length(miR100.val))),
  'mir_802_validated' = c(miR802.val, rep(NA, nr-length(miR802.val))),
  
  'mir_125a_predicted' = c(miR125a.pred, rep(NA, nr-length(miR125a.pred))),
  'mir_125b_predicted' = c(miR125b.pred, rep(NA, nr-length(miR125b.pred))),
  'mir_125_predicted' = c(miR125.pred, rep(NA, nr-length(miR125.pred))),
  'let_7c_predicted' = c(miRlet7c.pred, rep(NA, nr-length(miRlet7c.pred))),
  'let_7a_predicted' = c(miRlet7a.pred, rep(NA, nr-length(miRlet7a.pred))),
  'let_7e_predicted' = c(miRlet7e.pred, rep(NA, nr-length(miRlet7e.pred))),
  'mir_99a_predicted' = c(miR99a.pred, rep(NA, nr-length(miR99a.pred))),
  'mir_99b_predicted' = c(miR99b.pred, rep(NA, nr-length(miR99b.pred))),
  'mir_99_predicted' = c(miR99.pred, rep(NA, nr-length(miR99.pred))),
  'mir_155_predicted' = c(miR155.pred, rep(NA, nr-length(miR155.pred))),
  'mir_100_predicted' = c(miR100.pred, rep(NA, nr-length(miR100.pred))),
  'mir_802_predicted' = c(miR802.pred, rep(NA, nr-length(miR802.pred))),
  
  'mir_125_henning' = c(miR125.henning, rep(NA, nr-length(miR125.henning))),
  'let_7c_henning' = c(miRlet7c.henning, rep(NA, nr-length(miRlet7c.henning))),
  'mir_99_henning' = c(miR99.henning, rep(NA, nr-length(miR99.henning))),
  
  stringsAsFactors = FALSE
)
  
# add description row with na's
gmx <- rbind(rep("na", length(gmx)), gmx)

write.table(gmx, file="/mnt/projects/iamp/results/gsea/custom/mir_target_genesets.gmx", na = "", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
