# Jinliang Yang
# updated: 8/14/2014
# Purpose: zoomin plot

source("~/Documents/SAM_GS/SAM_proj/lib/getBayes.r")
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


######################
source("~/Documents/SAM_GS/SAM_proj/lib/findGenePhytozome.R")
nrow(subset(bayes3, ModelFreq > 0.02))
#21
sig <- subset(bayes3, ModelFreq > 0.02)
sig$chr <- paste("chr", sig$chr, sep="")

gene1 <- findGenePhytozome(fgs=NULL, sig=sig, binsize=100000, annotation=TRUE)
#FGSv2 loaded automatically!
#[ 21 ] SNPs > 0.02 cutoff
#[ 67 ] unique genes were identified with binsize [ 100000 ]
#[ 120 ] Phytozome annotations [Zmays_284_6a] were found!
write.table(gene1, "~/Documents/SAM_GS/SAM_proj/reports/samvolume_genes_100kbin_bayes.txt",
            sep="\t", row.names=FALSE, quote=FALSE)


