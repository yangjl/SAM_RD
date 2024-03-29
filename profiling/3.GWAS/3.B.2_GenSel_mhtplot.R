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

source("~/Documents/Rcodes/quickMHTplot.R")
par(mfrow=c(1,1))
#### pai=0.995
bayes1 <- getBayes("/home/NSF-SAM-GS/GenSel/SAM_run41000.mrkRes1")
#bayes1s <- subset(bayes1, ModelFreq > 0.005)
quickMHTplot(res=bayes1, main="SAM=380, pai=0.995", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")
  
#### pai=0.999
bayes2 <- getBayes("/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.mrkRes2")
#bayes2s <- subset(bayes2, ModelFreq > 0.005)
quickMHTplot(res=bayes2, main="SAM=380, pai=0.999, SNP~700", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")

#### pai=0.9995
bayes3 <- getBayes("/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.mrkRes1")
bayes3s <- subset(bayes3, ModelFreq > 0.01)
bayes3s <- bayes3s[order(bayes3s$ModelFreq, decreasing=TRUE),]

quickMHTplot(res=bayes3, main="SAM=380, pai=0.9995, SNP~350", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")

