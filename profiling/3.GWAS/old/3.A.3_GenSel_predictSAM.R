# Jinliang Yang
# Purpose: Prediction of SAM large and SAM small
# date: July.14.2014
# location: server.9

change_pheno <- function(){
  setwd("/home/NSF-SAM-GS/snpmatrix")
  pheno_large <- read.table("SAM_large_chrall.gensel.pheno", header=TRUE)
  names(pheno_large)[2] <- "log2vol"
  write.table(pheno_large[, -4], "/home/NSF-SAM-GS/GenSel/predictSAM/SAM_large.pheno", sep="\t",
              row.names=FALSE, quote=FALSE)
  
  pheno_small <- read.table("SAM_small_chrall.gensel.pheno", header=TRUE)
  names(pheno_small)[2] <- "log2vol"
  write.table(pheno_small[, -4], "/home/NSF-SAM-GS/GenSel/predictSAM/SAM_small.pheno", sep="\t",
              row.names=FALSE, quote=FALSE)
}
########
change_pheno()

### common files
source("~/Documents/Rcodes/sh4Predict.R")
pwd="/home/NSF-SAM-GS/GenSel/predictSAM/"
map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map"
mymrkRes <- "/home/NSF-SAM-GS/GenSel/SAM_run41000_pi.mrkRes1"

### prediction for large SAM
sh1 <- "SAM_predict_large.sh"
pheno1 <- "/home/NSF-SAM-GS/GenSel/predictSAM/SAM_large.pheno"
geno_large = "/home/NSF-SAM-GS/snpmatrix/SAM_large_chrall.gensel" 

sh4Predict(pwd=pwd, sh=sh1, geno=geno_large, pheno=pheno1, map=map, mrkRes=mymrkRes) 

### prediction for small SAM
sh2 <- "SAM_predict_small.sh"
pheno2 <- "/home/NSF-SAM-GS/GenSel/predictSAM/SAM_small.pheno"
geno_small = "/home/NSF-SAM-GS/snpmatrix/SAM_small_chrall.gensel" 

sh4Predict(pwd=pwd, sh=sh2, geno=geno_small, pheno=pheno2, map=map, mrkRes=mymrkRes) 




