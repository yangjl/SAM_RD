# Jinliang Yang
# updated: 7/14/2014
# Purpose: quick plot of GWAS results

getBayes <- function(inputfile="/home/NSF-SAM-GS/GenSel/SAM_run41000.mrkRes1"){
  res <- read.table(inputfile, header=TRUE)
  message(sprintf("input [ %s ]", nrow(res)))
  res$chr <- gsub("_.*", "", res$snpid)
  res$chr <- as.numeric(as.character(sub("chr", "", res$chr)))
  res <- subset(res, !is.na(chr))
  res$pos <- as.numeric(as.character(gsub(".*_", "", res$snpid)))
  message(sprintf("remove unknow chr, remainning [ %s ]", nrow(res)))
  return(res)
}

#### pai=0.9995
bayes3 <- getBayes("/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.mrkRes1")
#bayes2s <- subset(bayes2, ModelFreq > 0.005)

source("~/Documents/Rcodes/quickMHTplot.R")
quickMHTplot(res=bayes3, main="SAM=380, pai=0.9995, SNP~350", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")
abline(h=0.02, lwd=2, col="red", lty=2)

##############################################################################
subset(bayes3, ModelFreq > 0.02)
fgs <- read.table("~/bin/NGSbin/fgsv2.txt", header=TRUE);

subset(fgs, chr=="chr4" & start < 237392971 & end > 237392971)
#21744 GRMZM5G851807 chr4 237390273 237393310      +
subset(fgs, chr=="chr5" & start < 42451940-300000 & end > 42451940)
#23208 GRMZM2G014813 chr5 42530972 42536694      -
subset(fgs, chr=="chr5" & start < 208641840 & end > 208641840)
#26043 GRMZM2G380432 chr5 208640720 208642406      -


pval <- bayes3[, c("snpid", "chr", "pos", "ModelFreq")]
pval$chr <- paste("chr", pval$chr, sep="")
gene <- fgs

source("~/Documents/SAM_GS/SAM_proj/lib/Zoomin_plot.R")
ZoomInPlot(pval=pval, gene=gene, chr="chr4", window=20000, 
           snppos=237393027, ylim=c(0,0.06), idxmax=NULL)

ZoomInPlot(pval=pval, gene=gene, chr="chr5", window=20000, 
           snppos=42451940, ylim=c(0,0.10), idxmax=NULL)

ZoomInPlot(pval=pval, gene=gene, chr="chr5", window=10000, 
           snppos=208641840, ylim=c(0,0.04), idxmax=NULL)


