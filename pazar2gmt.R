p <- read.delim("/mnt/projects/generic/data/pazar/pazar_all.csv")
ph <- p[p$species=="Homo sapiens",]

library("biomaRt")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75

# translate TF ID into HGNC symbol
transcripts <- getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol"), mart=mart)
pha <- merge(ph, transcripts, by.x="Ensembl_TF_ID", by.y="ensembl_transcript_id")
names(pha)[names(pha)=="hgnc_symbol"] <- "TF_NAME"

# translate target gene ID into HGNC symbol
genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart=mart)
pha <- merge(pha, genes, by.x="Ensembl_Target_ID", by.y="ensembl_gene_id")
names(pha)[names(pha)=="hgnc_symbol"] <- "TARGET"

tf.melted <- unique(pha[,c("PAZAR_TF_ID", "TF_NAME", "TARGET")])
tf.melted <- tf.melted[!is.na(tf.melted$TF_NAME) & tf.melted$TF_NAME != "" & !is.na(tf.melted$TARGET) & tf.melted$TARGET != "",]

library(plyr)
tf.collapsed <- ddply(tf.melted, .(TF_NAME), summarize, TARGETS = paste(unique(TARGET), collapse = "\t"), DESC = paste(unique(PAZAR_TF_ID), collapse = ","))
tf.collapsed$TF_NAME <- paste0(tf.collapsed$TF_NAME, "_PAZAR")

write.table(tf.collapsed[c("TF_NAME", "DESC", "TARGETS")], file="/mnt/projects/generic/data/pazar/pazar.gmt", quote=F, sep="\t", row.names=F, col.names=F, na="")
