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
samGWASplot <- function(trait=traits[1], mmf = 0.001){
  setwd("/home/NSF-SAM-GS/GenSel/realrun/")
  source("~/Documents/Rcodes/quickMHTplot.R")  
  #### pai=0.995
  infile <- paste(trait, "_run41000.mrkRes1", sep="")
  bayes1 <- getBayes(infile)
  bayes2 <- subset(bayes1, ModelFreq > mmf)
  message(sprintf("after filtering [ %s ]", nrow(bayes2)))
  quickMHTplot(res=bayes2, main=trait, cex=.6, pch=19, 
               col=rep(c("slateblue", "cyan4"), 5), 
               GAP=5e+06, ylab="Model Frq", ylim=c(0, 0.08), 
               yaxis=c(0, 0.02, 0.04, 0.06, 0.08),
               col2plot="ModelFreq")
  abline(h=0.01, lty=2, col="red", lwd=2)
  return(bayes1)
}
#########################
traits <- c("Height", "Radius", "arc_length",
            "HtRa", "area", "volume",  "para_coeff", 
             "surface_area")

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


##################
bayes1$trait <- traits[1]
bayes2$trait <- traits[2]
bayes3$trait <- traits[3]
bayes4$trait <- traits[4]
bayes5$trait <- traits[5]
bayes6$trait <- traits[6]
bayes7$trait <- traits[7]
bayes8$trait <- traits[8]

bayesall <- rbind(bayes1, bayes2, bayes3, bayes4, bayes5, bayes6, bayes7, bayes8)
#########################
bayes0 <- bayesall[, c("trait", "snpid", "chr", "pos", "ModelFreq","Effect" )]
tav <- subset(bayes0, ModelFreq > 0.01)
table(tav$trait)

write.table(bayes0, "~/Documents/SAM_GS/SAM_proj/reports/Bayes_8traits_allsnps.csv",
            sep=",", row.names=FALSE, quote=FALSE)
write.table(tav, "~/Documents/SAM_GS/SAM_proj/reports/Bayes_8traits_TAVs.csv",
            sep=",", row.names=FALSE, quote=FALSE)
