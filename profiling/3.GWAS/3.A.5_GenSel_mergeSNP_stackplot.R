# Jinliang Yang
# updated: 8/14/2014
# Purpose: quick plot of GWAS results

###########################
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
###########################
samGWASplot <- function(trait=traits[1], mmf = 0.001, plotlayout=TRUE){
  setwd("/home/NSF-SAM-GS/GenSel/realrun/")
  
  #### pai=0.995
  infile <- paste(trait, "_run41000.mrkRes1", sep="")
  bayes1 <- getBayes(infile)
  bayes2 <- subset(bayes1, ModelFreq > mmf)
  message(sprintf("after filtering [ %s ]", nrow(bayes2)))
  source("~/Documents/SAM_GS/SAM_proj/lib/quickMHTplot_stacking.R")  
  if(plotlayout){
    quickMHTplot(res=bayes2, main=trait, cex=.6, pch=19, 
                 col=rep(c("slateblue", "cyan4"), 5), 
                 GAP=5e+06, ylab="Model Frq", ylim=c(0, 0.08), 
                 yaxis=c(0, 0.02, 0.04, 0.06, 0.08),
                 col2plot="ModelFreq")
  }
  
  for(i in 1:10){
    points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
           pch = pch, col=col[i], cex=cex);
  }
  
  abline(h=0.01, lty=2, col="red", lwd=2)
  return(bayes1)
}
#########################
traits <- c("arc_length", "area", "Height", "HtRa", "para_coeff", 
            "Radius", "surface_area", "volume")

pdf("~/Documents/SAM_GS/SAM_proj/graphs/Bayes_plot.pdf", width=8, height=12)
par(mfrow=c(8,1), mar=c(3,4,2,2), oma=c(0, 3, 0, 3))
bayes1 <- samGWASplot(trait=traits[1], mmf = 0.001)  
bayes2 <- samGWASplot(trait=traits[2], mmf = 0.001)  
bayes3 <- samGWASplot(trait=traits[3], mmf = 0.001)  
bayes4 <- samGWASplot(trait=traits[4], mmf = 0.001)  
bayes5 <- samGWASplot(trait=traits[5], mmf = 0.0008)  
bayes6 <- samGWASplot(trait=traits[6], mmf = 0.001)  
bayes7 <- samGWASplot(trait=traits[7], mmf = 0.001)  
bayes8 <- samGWASplot(trait=traits[8], mmf = 0.001)  
dev.off()





















#########################
bayes3 <- bayes3[order(bayes3$ModelFreq, decreasing=TRUE),]
bayes3 <- bayes3[,c("snpid", "chr", "pos", "ModelFreq", "Effect")]

write.table(bayes3, "~/Documents/SAM_GS/SAM_proj/reports/bayesC_860K_pai9995.csv",
            sep=",", row.names=FALSE, quote=FALSE)
