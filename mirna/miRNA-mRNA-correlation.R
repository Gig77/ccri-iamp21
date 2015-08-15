mirna <- read.delim("/mnt/projects/iamp/data/henning/mirna-expression-henning.tsv", check.names = F, stringsAsFactors = F)

# set corresponding CCRI sample IDs

mirna$ccri <- NA
mirna$ccri[mirna$'Sample #' == "3"] <- "C5"
mirna$ccri[mirna$'Sample #' == "4"] <- "C6"
mirna$ccri[mirna$'Sample #' == "10"] <- "C1"
mirna$ccri[mirna$'Sample #' == "12"] <- "C2"
mirna$ccri[mirna$'Sample #' == "13"] <- "C3"
mirna$ccri[mirna$'Sample #' == "14"] <- "C4"
mirna$ccri[mirna$'Sample #' == "17"] <- "C7"
mirna$ccri[mirna$'Sample #' == "19"] <- "C8"
mirna$ccri[mirna$'Sample #' == "23"] <- "C9"
mirna$ccri[mirna$'Sample #' == "25"] <- "C10"
mirna$ccri[mirna$'Sample #' == "27"] <- "C11"
mirna$ccri[mirna$'Sample #' == "29"] <- "C12"
mirna$ccri[mirna$'Sample #' == "30"] <- "A1"
mirna$ccri[mirna$'Sample #' == "31"] <- "A2"
mirna$ccri[mirna$'Sample #' == "32"] <- "A3"
mirna$ccri[mirna$'Sample #' == "33"] <- "A4"
mirna$ccri[mirna$'Sample #' == "34"] <- "A5"
mirna$ccri[mirna$'Sample #' == "35"] <- "A6"
mirna$ccri[mirna$'Sample #' == "36"] <- "A7"
mirna$ccri[mirna$'Sample #' == "37"] <- "A8"
mirna$ccri[mirna$'Sample #' == "38"] <- "A9"
mirna$ccri[mirna$'Sample #' == "39"] <- "A10"
mirna$ccri[mirna$'Sample #' == "40"] <- "A11"
mirna$ccri[mirna$'Sample #' == "41"] <- "A12"
mirna$ccri[mirna$'Sample #' == "42"] <- "A13"
mirna$ccri[mirna$'Sample #' == "D1"] <- "D1"
mirna$ccri[mirna$'Sample #' == "D2"] <- "D2"
mirna$ccri[mirna$'Sample #' == "D3"] <- "D3"
mirna$ccri[mirna$'Sample #' == "D4"] <- "D4"
mirna$ccri[mirna$'Sample #' == "D5"] <- "D5"
mirna$ccri[mirna$'Sample #' == "D6"] <- "D6"
mirna$ccri[mirna$'Sample #' == "D7"] <- "D7"
mirna$ccri[mirna$'Sample #' == "D8"] <- "D8"
mirna$ccri[mirna$'Sample #' == "S1"] <- "S1"
mirna$ccri[mirna$'Sample #' == "S2"] <- "S2"
mirna$ccri[mirna$'Sample #' == "S3"] <- "S3"

# factor for leukemia subtype

mirna$disease <- NA
mirna$disease[grepl("^A", mirna$ccri)] <- "iAMP21"
mirna$disease[grepl("^C", mirna$ccri)] <- "E/R"
mirna$disease[grepl("^D", mirna$ccri)] <- "DS"
mirna$disease[grepl("^S", mirna$ccri)] <- "B-cell"
mirna$disease <- factor(mirna$disease, levels = c("E/R", "iAMP21", "DS", "B-cell"))

# remove failed samples

mirna <- mirna[!mirna$ccri %in% c("C2", "A4"),]

# get normalized mRNA gene expression

library("DESeq2")
countFiles <- read.delim("/mnt/projects/iamp/results/anduril/execute/_deseqExprMatrix_geneCounts_array1/array/_index")
cds <- DESeqDataSetFromHTSeqCount(sampleTable=countFiles, directory="/", design=~1)
expressed <- rowSums(counts(cds)) >= 20 | apply(counts(cds), 1, max) >= 10
invariant <- apply(counts(cds), 1, sd) < 3
vsd <- varianceStabilizingTransformation(cds[expressed & !invariant,])
mat <- assay(vsd)

# translate into HGNC gene names (if known)
library("biomaRt")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart=mart)
hgnc <- aggregate(hgnc_symbol~ensembl_gene_id, paste, collapse=",", data=hgnc)
id2name <- merge(data.frame(ensembl_gene_id=rownames(mat), stringsAsFactors = F), hgnc, all.x=T)
id2name$hgnc_symbol[id2name$hgnc_symbol == ""] <- id2name$ensembl_gene_id[id2name$hgnc_symbol == ""]
rownames(mat) <- id2name$hgnc_symbol

# remove samples without measurement
mat <- mat[,!colnames(mat) %in% c("C2", "A4")]

# add miR expression to matrix
mirmat <- mirna[!is.na(mirna$ccri), c("miR-99a dCt", "let-7c dCt", "miR-125b dCt", "miR-100 dCt", "let-7a dCt", "miR-99b dCt", "let-7e dCt", "miR-125a dCt",
                                      "miR-155 dCt", "miR-802 dCt")]
colnames(mirmat) <- sub(" dCt", "", colnames(mirmat))
rownames(mirmat) <- mirna$ccri[!is.na(mirna$ccri)]
mirmat <- t(mirmat)
mirmat <- mirmat[,colnames(mat)]
mirmat <- log(mirmat) # we look for anti-correlated genes (miR down -> target up)

# compute correlation
corr <- cor(t(mat), t(mirmat), method="spearman")

# read validated/predicted miR target genes
validated <- read.delim("/mnt/projects/iamp/data/miRecords/miRecords.validatedTargets.txt", stringsAsFactors = FALSE)
predicted <- read.delim("/mnt/projects/iamp/data/miRecords/miRecords.predictedTargets.txt", stringsAsFactors = FALSE)


# plot results
sample.colors = rep(NA, length(colnames(mat)))
sample.colors[grepl("^S", colnames(mat))] <- "darkgray"
sample.colors[grepl("^A", colnames(mat))] <- "blue"
sample.colors[grepl("^D", colnames(mat))] <- "red"
sample.colors[grepl("^C", colnames(mat))] <- "black"

genesPerPage <- 16
R.cutoff <- 0.7
p.cutoff <- 1e-5

pdf("/mnt/projects/iamp/results/mirna/miR-mRNA-top-correlated-genes.pdf", height=10, width=10)

# plot pairwise correlations among miRs

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}
pairs(t(mirmat),lower.panel=panel.smooth, upper.panel=panel.cor)

# plot correlations miRs <--> other genes

for(m in c(1:length(colnames(corr)))) {
  mir.name <- colnames(corr)[m]
  for (type in c("correlated", "anti-correlated")) {
    if (type == "correlated") {
      selected <- rownames(corr)[corr[,m] >= R.cutoff]
      selected <- selected[order(corr[selected,m], decreasing = T)]
    } else {
      selected <- rownames(corr)[corr[,m] <= -R.cutoff]
      selected <- selected[order(corr[selected,m], decreasing = F)]
    }

    # throw out non-significant correlations before plotting
    selected.keep <- character(0)
    for(s in selected) {
      test <- cor.test(mirmat[mir.name,],mat[s,],method="spearman")
      if (test$p.value <= p.cutoff) selected.keep <- c(selected.keep, s)
    }
    
    # no significant associations left?
    if (length(selected.keep) == 0) next
    
    # plot page by page
    for (page in 1:ceiling(length(selected.keep)/genesPerPage)) {
      par(mfrow=c(4,4), mar=c(2,2,2,1), oma=c(0,0,2,0))
      for(i in ((page-1)*genesPerPage+1):min(page*genesPerPage,length(selected.keep))) {
        gene.name <- selected.keep[i]
        x <- mirmat[mir.name,]
        y <- mat[gene.name,]
        fit <- lm(y~x)
        test <- cor.test(x,y,method="spearman")
        plot(x, y, xlab = "", ylab = "", ylim=c(min(0, y), max(y)), xlim=c(min(0, x), max(x)), pch = "")
        
        # highlight validated or predicted target genes
        suffix <- ""
        if (sum(grepl(mir.name, validated$miRNA_mature_ID) & gene.name == validated$symbol) > 0) {
          print(paste("VALIDATED", mir.name, gene.name))
          suffix <- "**"
        } else if (sum(grepl(mir.name, predicted$miRNA.ID) & gene.name == predicted$Symbol) > 0) {
          print(paste("PREDICTED", mir.name, gene.name))
          suffix <- "*"
        }
        mtext(text = paste0(gene.name, suffix), line = 0, cex = 0.8)
        text(x, y, colnames(mat), col = sample.colors, cex = 0.8)
        abline(fit, col="red")
        text((min(x)+max(x))/2, min(0, y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)
      }
      mtext(paste0("Genes ", type, " in expression with ", mir.name, " (|R|>=", R.cutoff, ", p<=", p.cutoff, ") [page ", page, " of ", ceiling(length(selected.keep)/genesPerPage), "]"), outer=T)
    }
  }
}
dev.off()

           