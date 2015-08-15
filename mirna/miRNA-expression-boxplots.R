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

# remove outlier

mirna <- mirna[!mirna$ccri %in% c("C2", "A4"),]

doplot <- function(mir) {
  expr <- mirna[,paste(mir, "dCt")] * 100
  boxplot(expr~mirna$disease, log="y", ylim=c(1e-4, 100), yaxt="n", ylab=paste(mir, "(% of U6)"), main=mir, na.action=na.exclude, outline=F)
  stripchart(expr~mirna$disease, log="y", ylim=c(1e-4, 100), yaxt="n", method="jitter", jitter=0.3, na.action=na.exclude, vertical=T, pch=19, col=1:length(levels(as.factor(mirna$disease))), add=T)
  abline(h=c(1e-4, 1e-3, 0.01, 0.1, 1, 10, 100), lty=2, col="darkgray")
  axis(2, at=c(1e-4, 1e-3, 0.01, 0.1, 1, 10, 100), label=c("1e-4", "1e-3", "0.01", "0.1", "1", "10", "100"), las=2)
  axis(2, at=c(seq(0.0002, 0.0009, 0.0001), seq(0.002, 0.009, 0.001), seq(0.02, 0.09, 0.01), seq(0.2, 0.9, 0.1), seq(2, 9, 1), seq(20, 90, 10)), label=F)
}

pdf("/mnt/projects/iamp/results/mirna/mirna-expression-boxplots.pdf", height=12)
par(mfrow=c(4,3), mar=c(3,4,1,1))
doplot("miR-155")
doplot("miR-99a")
doplot("miR-100")
doplot("miR-99b")
doplot("miR-802")
doplot("let-7c")
doplot("let-7a")
doplot("let-7e")
doplot("miR-125b")
doplot("miR-125a")
dev.off()


