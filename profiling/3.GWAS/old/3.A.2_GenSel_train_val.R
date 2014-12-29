# Jinliang Yang
# Purpose: GenSel of cross-validation 80% training and 20% val
# date: July.14.2014
# location: server.9


pwd="/home/NSF-SAM-GS/GenSel/train_val"
geno="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.gensel.newbin" 
map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map"

source("~/Documents/Rcodes/sh4GenSel.R")
for(i in 1:10){
  mysh <- paste("SAM_run41000_s", i, ".sh", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_s", i, "_train.txt", sep="")
  sh4GenSel(pwd=pwd, sh=mysh, pi=0.9995, geno=geno, map=map, pheno=mypheno,        
            chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)  
}


### prediction
source("~/Documents/Rcodes/sh4Predict.R")
for(i in 1:10){
  mysh <- paste("SAM_predict_s", i, ".sh", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_s", i, "_val.txt", sep="")
  mymrkRes <- paste("/home/NSF-SAM-GS/GenSel/train_val/SAM_run41000_s", i, ".mrkRes1", sep="")
  sh4Predict(pwd=pwd, sh=mysh, geno=geno, pheno=mypheno, map=map, mrkRes=mymrkRes) 
                        
}


