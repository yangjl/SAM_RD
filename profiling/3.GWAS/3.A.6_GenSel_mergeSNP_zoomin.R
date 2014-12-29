# Jinliang Yang
# updated: 8/14/2014
# Purpose: zoomin plot

#### pai=0.9995
bayes3 <- getBayes("/media/RAID/NSF-SAM-GS/GenSel/run_mergedSNPs/SAM_run41000_pi9995.mrkRes1")
#bayes3s <- subset(bayes3, ModelFreq > 0.01)
#bayes3s <- bayes3s[order(bayes3s$ModelFreq, decreasing=TRUE),]

quickMHTplot(res=bayes3, main="SAM=380, pai=0.9995, SNP~350", cex=.5, pch=19, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=1e+07, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")
abline(h=0.02, lty=2, col="red", lwd=3)


bayes3 <- bayes3[order(bayes3$ModelFreq, decreasing=TRUE),]

##############################################################################
nrow(subset(bayes3, ModelFreq > 0.02))
#21

sig <- subset(bayes3, ModelFreq > 0.02)
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

##############
source("~/Documents/SAM_GS/SAM_proj/lib/Zoomin_plot.R")
chr2 <- ZoomInPlot(pval=pval, gene=gene, chr="chr2", window=500000, 
                   snppos=79480554, ylim=c(0,0.08), idxmax=NULL)
abline=c(h=0.02, col="red", lty=2, lwd=2)

chr4 <- ZoomInPlot(pval=pval, gene=gene, chr="chr4", window=10000, 
                   snppos=237392971, ylim=c(0,0.08), idxmax=NULL)
abline=c(h=0.02, col="red", lty=2, lwd=2)

chr5 <- ZoomInPlot(pval=pval, gene=gene, chr="chr5", window=10000, 
                   snppos=72384004, ylim=c(0,0.09), idxmax=NULL)
abline=c(h=0.02, col="red", lty=2, lwd=2)

chr8 <- ZoomInPlot(pval=pval, gene=gene, chr="chr8", window=10000, 
                   snppos=72384004, ylim=c(0,0.09), idxmax=NULL)
abline=c(h=0.02, col="red", lty=2, lwd=2)

