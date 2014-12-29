# Jinliang Yang
# updated: 8/14/2014
# Purpose: quick plot of GWAS results


###########################
samWINplot <- function(trait=traits[1]){
  setwd("/home/NSF-SAM-GS/GenSel/realrun/")
  source("~/Documents/Rcodes/quickMHTplot.R")
  par(mfrow=c(1,1))
  #### pai=0.995
  infile <- paste(trait, "_run41000.winQTL1", sep="")
  bayes1 <- read.table(infile, header=FALSE, skip=1)
  names(bayes1) <- c("window", "start", "end", "nosnps", "var", "cumvar", "p", "pavg",
                     "mappos0", "mapposn", "chr_Mb")
  bayes1$chr <- gsub("_.*", "", bayes1$chr_Mb)
  bayes1$pos <- as.numeric(as.character(gsub(".*_", "", bayes1$chr_Mb)))*1000000
  quickMHTplot(res=bayes1, main=trait, cex=.6, pch=19, 
               col=rep(c("slateblue", "cyan4"), 5), 
               GAP=5e+06, ylab="model frequency", yaxis=NULL,
               col2plot="var")
  return(bayes1)
}
#########################
traits <- c("arc_length", "area", "Height", "HtRa", "para_coeff", 
            "Radius", "surface_area", "volume")
bayes1 <- samWINplot(trait=traits[1])  
bayes2 <- samWINplot(trait=traits[2])  
bayes3 <- samWINplot(trait=traits[3])  
bayes4 <- samWINplot(trait=traits[4])  
bayes5 <- samWINplot(trait=traits[5])  
bayes6 <- samWINplot(trait=traits[6])  
bayes7 <- samWINplot(trait=traits[7])  
bayes8 <- samWINplot(trait=traits[8])  























#########################
bayes3 <- bayes3[order(bayes3$ModelFreq, decreasing=TRUE),]
bayes3 <- bayes3[,c("snpid", "chr", "pos", "ModelFreq", "Effect")]

write.table(bayes3, "~/Documents/SAM_GS/SAM_proj/reports/bayesC_860K_pai9995.csv",
            sep=",", row.names=FALSE, quote=FALSE)
