s <- read.delim("/mnt/projects/iamp/results/qc/allpatients.stats.txt")

rownames(s) <- sub(".*C57C3ACXX_CV_([^_]+)_.*", "\\1", s$sample)

s$type <- NA
s$type[grepl("_C\\d", s$sample)] <- "ER"
s$type[grepl("_A\\d", s$sample)] <- "iAMP21"
s$type[grepl("_D\\d", s$sample)] <- "DS"
s$type[grepl("_S\\d", s$sample)] <- "CD19+"
s$type <- as.factor(s$type)

s$pct_rrna <- (s$exonic-s$non.rRNA)/s$exonic
s$pct_rna <- (s$non.rRNA-s$protein)/s$non.rRNA
s$pct_intronic <- (s$mapped-s$exonic)/s$mapped
s$pct_dup <- (s$uniquely.mapped-s$non.duplicates)/s$uniquely.mapped

test <- t.test(s$pct_rrna[s$type %in% c("ER", "iAMP21")], s$pct_rrna[s$type %in% c("DS")])

pdf("/mnt/projects/iamp/results/pct-rrna-per-subtype.pdf")
boxplot(pct_rrna~type, data=s, xlab=sprintf("Subtype\np(DS vs. ER+iAMP21ER) = %.2g", test$p.value), ylab="Percent rRNA in total exonic reads", main="Percentage of rRNA-reads per subtype", outline=F, cex.axis=0.8)
stripchart(s$pct_rrna~s$type, method="jitter", vertical=T, pch=19, col=1:length(levels(as.factor(s$type))), add=T)
dev.off()

s[order(s$pct_rrna, decreasing=TRUE), c("sample", "pct_rrna")]

plot(hclust(dist(s[,c("pct_rrna", "pct_rna", "pct_intronic", "pct_dup")])))

pca <- prcomp(s[,c("pct_rrna", "pct_rna", "pct_intronic", "pct_dup")], center = TRUE, scale. = TRUE) 
par(mar=c(4,4,4,7))
par(xpd=TRUE)
plot( pca$x[,1], pca$x[,2], col="white", xlab="PC1", ylab="PC2", main="Sample clustering by library properties")
text( pca$x[,1], pca$x[,2], rownames(s), col=as.numeric(s$type) )
legend(x=par("usr")[2], y=max(par("usr")[3:4]), levels(s$type), col=1:length(levels(s$type)), pch=1)
