# read data
library("gdata")
cluster1 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_01", header = TRUE)
cluster2 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_02", header = TRUE)
cluster3 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_03", header = TRUE)
cluster4 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_04", header = TRUE)
cluster5 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_05", header = TRUE)
cluster6 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_06", header = TRUE)
cluster7 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_07", header = TRUE)
cluster8 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_08", header = TRUE)
cluster9 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_09", header = TRUE)
cluster10 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_10", header = TRUE)
cluster11 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_11", header = TRUE)
cluster12 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_12", header = TRUE)
cluster13 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_13", header = TRUE)
cluster14 = read.xls("/mnt/projects/iamp/data/dick_2013/table-S2.xls", sheet = "Kmeans_14_cluster_14", header = TRUE)

# assign clusters to cell types according to annotation (sheet 1 in above excel file)
cluster.mlp <- rbind(cluster1, cluster2, cluster3, cluster4)
cluster.prog <- rbind(cluster5, cluster6, cluster8, cluster10)
cluster.eryth <- cluster7
cluster.myel_lymph <- cluster9
cluster.lymph <- rbind(cluster11, cluster12, cluster14)
cluster.myel <- cluster13

# keep only top-n most highly expressed genes per cluster
cluster.mlp <- cluster.mlp[order(cluster.mlp$HSC, decreasing = T),][1:500,]
cluster.prog <- cluster.prog[order(rowMeans(cluster.prog[,c("HSC", "PROB")])),][1:500,]
cluster.eryth <- cluster.eryth[order(cluster.eryth$MEP, decreasing = T),][1:500,]
cluster.myel_lymph <- cluster.myel_lymph[order(rowMeans(cluster.myel_lymph[,c("HSC", "MLP", "GMP", "PROB")]), decreasing = T),][1:500,]
cluster.lymph <- cluster.lymph[order(cluster.lymph$PROB, decreasing = T),][1:500,]
cluster.myel <- cluster.myel[order(cluster.myel$GMP, decreasing = T),][1:500,]

# map Illumina probe IDs to HGNC gene symbols
library("biomaRt")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "illumina_humanht_12_v4"), mart=mart)
hgnc <- hgnc[hgnc$hgnc_symbol != "" & hgnc$illumina_humanht_12_v4 != "",]
hgnc <- aggregate(hgnc_symbol~illumina_humanht_12_v4, paste, collapse=",", data=hgnc)

cluster.mlp.hgnc <- merge(cluster.mlp, hgnc, by.x = "IlluminaProbe", by.y = "illumina_humanht_12_v4", all.x = T)
cluster.prog.hgnc <- merge(cluster.prog, hgnc, by.x = "IlluminaProbe", by.y = "illumina_humanht_12_v4", all.x = T)
cluster.eryth.hgnc <- merge(cluster.eryth, hgnc, by.x = "IlluminaProbe", by.y = "illumina_humanht_12_v4", all.x = T)
cluster.myel_lymph.hgnc <- merge(cluster.myel_lymph, hgnc, by.x = "IlluminaProbe", by.y = "illumina_humanht_12_v4", all.x = T)
cluster.lymph.hgnc <- merge(cluster.lymph, hgnc, by.x = "IlluminaProbe", by.y = "illumina_humanht_12_v4", all.x = T)
cluster.myel.hgnc <- merge(cluster.myel, hgnc, by.x = "IlluminaProbe", by.y = "illumina_humanht_12_v4", all.x = T)

# convert to .gmt format
gs.mlp <- paste(sort(unique(unlist(strsplit(paste(unique(cluster.mlp.hgnc$hgnc_symbol[!is.na(cluster.mlp.hgnc$hgnc_symbol)]), collapse = ","), ",")))), collapse = "\t")
gs.mlp.df <- data.frame(name="DICK_2013_STEM_MLP", desc="Dick (2013) MLP-specific", genes=gs.mlp)

gs.prog <- paste(sort(unique(unlist(strsplit(paste(unique(cluster.prog.hgnc$hgnc_symbol[!is.na(cluster.prog.hgnc$hgnc_symbol)]), collapse = ","), ",")))), collapse = "\t")
gs.prog.df <- data.frame(name="DICK_2013_PROGENITOR", desc="Dick (2013) Progenitor", genes=gs.prog)

gs.eryth <- paste(sort(unique(unlist(strsplit(paste(unique(cluster.eryth.hgnc$hgnc_symbol[!is.na(cluster.eryth.hgnc$hgnc_symbol)]), collapse = ","), ",")))), collapse = "\t")
gs.eryth.df <- data.frame(name="DICK_2013_ERYTHROID", desc="Dick (2013) Erythroid", genes=gs.eryth)

gs.myel_lymph <- paste(sort(unique(unlist(strsplit(paste(unique(cluster.myel_lymph.hgnc$hgnc_symbol[!is.na(cluster.myel_lymph.hgnc$hgnc_symbol)]), collapse = ","), ",")))), collapse = "\t")
gs.myel_lymph.df <- data.frame(name="DICK_2013_MYELOID_LYMPHOID", desc="Dick (2013) Myeloid-lymphoid", genes=gs.myel_lymph)

gs.lymph <- paste(sort(unique(unlist(strsplit(paste(unique(cluster.lymph.hgnc$hgnc_symbol[!is.na(cluster.lymph.hgnc$hgnc_symbol)]), collapse = ","), ",")))), collapse = "\t")
gs.lymph.df <- data.frame(name="DICK_2013_LYMPHOID", desc="Dick (2013) Lymphoid", genes=gs.lymph)

gs.myel <- paste(sort(unique(unlist(strsplit(paste(unique(cluster.myel.hgnc$hgnc_symbol[!is.na(cluster.myel.hgnc$hgnc_symbol)]), collapse = ","), ",")))), collapse = "\t")
gs.myel.df <- data.frame(name="DICK_2013_MYELOID", desc="Dick (2013) Myeloid", genes=gs.myel)

# write gmt output
gs.combined <- rbind(gs.mlp.df, gs.prog.df, gs.eryth.df, gs.myel_lymph.df, gs.lymph.df, gs.myel.df)
write.table(gs.combined, file="/mnt/projects/generic/data/laurenti_2013_hematopoietic_lineages.gmt", quote=F, sep="\t", row.names=F, col.names=F, na="")
