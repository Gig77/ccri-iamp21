library(reshape)

# read expression matrix
expr <- read.delim("/mnt/projects/iamp/results/anduril/execute/deseqExprMatrix/expr.csv", stringsAsFactors = F, row.names = 1)

#ds.top100up <-gsub("\"", "", unlist(strsplit(system('cat /mnt/projects/iamp/results/anduril/execute/deseq_DSvsNonDS/results.csv | grep -vP "\tNA\t" | grep -P "^[^\t]+\t[^-]" | sort -k 4g | cut -f 1 2>/dev/null | grep -v "\"ids\"" 2>/dev/null | head -100', intern=TRUE), " ")))
#ds.top100dn <-gsub("\"", "", unlist(strsplit(system('cat /mnt/projects/iamp/results/anduril/execute/deseq_DSvsNonDS/results.csv | grep -vP "\tNA\t" | grep -P "^[^\t]+\t-" | sort -k 4g | cut -f 1 2>/dev/null | grep -v "\"ids\"" 2>/dev/null | head -100', intern=TRUE), " ")))

ds.top100up <-gsub("\"", "", unlist(strsplit(system("cut -f 1 /mnt/projects/iamp/results/anduril/execute/degTableUp_DSvsNonDS-degsFiltered/csv.csv | grep -v Ensembl | head -100", intern=TRUE), " ")))
ds.top100dn <-gsub("\"", "", unlist(strsplit(system("cut -f 1 /mnt/projects/iamp/results/anduril/execute/degTableDn_DSvsNonDS-degsFiltered/csv.csv | grep -v Ensembl | head -100", intern=TRUE), " ")))

# add gene symbols

# compare GC content
#library(biomaRt)
#mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
#genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "percentage_gc_content"), mart=mart)
#up <- merge(ds.top100up, genes, by.x=1, by.y="ensembl_gene_id")
#up$group <- "up"
#dn <- merge(ds.top100dn, genes, by.x=1, by.y="ensembl_gene_id")
#dn$group <- "dn"
#updn <- rbind(up, dn)
#updn$group <- as.factor(updn$group)
#boxplot(percentage_gc_content~group, updn)

expr.top100up <- expr[rownames(expr) %in% ds.top100up,]
expr.top100dn <- expr[rownames(expr) %in% ds.top100dn,]

# ---- top 100 DS-upregulated genes in iAMP21 samples
expr.top100up.iamp <- expr.top100up[,c("A2", "A4", "A7", "A10", "A11", "A12", "A13", "A1", "A3", "A8", "A9")]
expr.top100up.iamp.scaled <- t(scale(t(expr.top100up.iamp)))

molten.top100up.iamp <- melt(expr.top100up.iamp.scaled)
names(molten.top100up.iamp)[1:2] <- c("hgnc", "sample")

molten.top100up.iamp$group <- NA
molten.top100up.iamp$group[molten.top100up.iamp$sample %in% c("A2", "A4", "A7")] <- "Conny (3 samples)"
molten.top100up.iamp$group[molten.top100up.iamp$sample %in% c("A10", "A11", "A12", "A13", "A1", "A3", "A8", "A9")] <- "Doris (8 samples)"
molten.top100up.iamp$group <- as.factor(molten.top100up.iamp$group)

# ---- top 100 DS-downregulated genes in iAMP21 samples
expr.top100dn.iamp <- expr.top100dn[,c("A2", "A4", "A7", "A10", "A11", "A12", "A13", "A1", "A3", "A8", "A9")]
expr.top100dn.iamp.scaled <- t(scale(t(expr.top100dn.iamp)))

molten.top100dn.iamp <- melt(expr.top100dn.iamp.scaled)
names(molten.top100dn.iamp)[1:2] <- c("hgnc", "sample")

molten.top100dn.iamp$group <- NA
molten.top100dn.iamp$group[molten.top100dn.iamp$sample %in% c("A2", "A4", "A7")] <- "Conny (3 samples)"
molten.top100dn.iamp$group[molten.top100dn.iamp$sample %in% c("A10", "A11", "A12", "A13", "A1", "A3", "A8", "A9")] <- "Doris (8 samples)"
molten.top100dn.iamp$group <- as.factor(molten.top100dn.iamp$group)

# ---- top 100 DS-upregulated genes in E/R samples
expr.top100up.ER<- expr.top100up[,c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")]
expr.top100up.ER.scaled <- t(scale(t(expr.top100up.ER)))

molten.top100up.ER <- melt(expr.top100up.ER.scaled)
names(molten.top100up.ER)[1:2] <- c("hgnc", "sample")

molten.top100up.ER$group <- NA
molten.top100up.ER$group[molten.top100up.ER$sample %in% c("C4", "C6")] <- "Conny (2 samples)"
molten.top100up.ER$group[molten.top100up.ER$sample %in% c("C1", "C2", "C3", "C5", "C7", "C8", "C9", "C10", "C11", "C12")] <- "Doris (10 samples)"
molten.top100up.ER$group <- as.factor(molten.top100up.ER$group)

# ---- top 100 DS-downregulated genes in E/R samples
expr.top100dn.ER <- expr.top100dn[,c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")]
expr.top100dn.ER.scaled <- t(scale(t(expr.top100dn.ER)))

molten.top100dn.ER <- melt(expr.top100dn.ER.scaled)
names(molten.top100dn.ER)[1:2] <- c("hgnc", "sample")

molten.top100dn.ER$group <- NA
molten.top100dn.ER$group[molten.top100dn.ER$sample %in% c("C4", "C6")] <- "Conny (2 samples)"
molten.top100dn.ER$group[molten.top100dn.ER$sample %in% c("C1", "C2", "C3", "C5", "C7", "C8", "C9", "C10", "C11", "C12")] <- "Doris (10 samples)"
molten.top100dn.ER$group <- as.factor(molten.top100dn.ER$group)

# --- plot side-by-side

test.up.iamp <- kruskal.test(molten.top100up.iamp$value, molten.top100up.iamp$group, na.action=na.exclude)
test.up.ER <- kruskal.test(molten.top100up.ER$value, molten.top100up.ER$group, na.action=na.exclude)
test.dn.iamp <- kruskal.test(molten.top100dn.iamp$value, molten.top100dn.iamp$group, na.action=na.exclude)
test.dn.ER <- kruskal.test(molten.top100dn.ER$value, molten.top100dn.ER$group, na.action=na.exclude)

pdf("/mnt/projects/iamp/results/expr-DS-DEGs-in-iAMP21-and-ER-samples.pdf", width=12, height=12)
par(mfrow=c(2,2))

boxplot(value~group, molten.top100up.iamp, xlab="RNA extracted from", ylab="Normalized expression (Z-score)", main="Expression of top-100 DS-upregulated genes in iAMP21 samples")  
text(1.5, 2.3, sprintf("p(Kruskal) = %.2g", test.up.iamp$p.value))

boxplot(value~group, molten.top100up.ER, xlab="RNA extracted from", ylab="Normalized expression (Z-score)", main="Expression of top-100 DS-upregulated genes in E/R samples")  
text(1.5, 2, sprintf("p(Kruskal) = %.2g", test.up.ER$p.value))

boxplot(value~group, molten.top100dn.iamp, xlab="RNA extracted from", ylab="Normalized expression (Z-score)", main="Expression of top-100 DS-downregulated genes in iAMP21 samples")  
text(1.5, 1.5, sprintf("p(Kruskal) = %.2g", test.dn.iamp$p.value))

boxplot(value~group, molten.top100dn.ER, xlab="RNA extracted from", ylab="Normalized expression (Z-score)", main="Expression of top-100 DS-downregulated genes in E/R samples")  
text(1.5, 1.5, sprintf("p(Kruskal) = %.2g", test.dn.ER$p.value))

dev.off()

