# Jinliang Yang
# Purpose: GenSel of cross-validation 80% training and 20% val
# date: July.14.2014
# location: server.9


pwd="/home/NSF-SAM-GS/GenSel/rc_runs"
geno="/home/NSF-SAM-GS/GenSel/rc_runs/SAM_trc_rpkm_rescaled.txt.newbin" 
map="/home/NSF-SAM-GS/GenSel/rc_runs/SAM_map.txt"

source("~/Documents/Rcodes/sh4GenSel.R")
### test run
mysh <- paste("SAMrc_test11000_s", i, ".sh", sep="")
mypheno <- paste("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_s", 2, "_train.txt", sep="")
sh4GenSel(pwd=pwd, sh=mysh, pi=0.995, geno=geno, map=map, pheno=mypheno, findsale ="no",       
          chainLength=11000, burnin=1000, varGenotypic=1, varResidual=1)  

######### run41000
for(i in 1:10){
  mysh <- paste("SAMrc_run41000_s", i, ".sh", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_s", i, "_train.txt", sep="")
  sh4GenSel(pwd=pwd, sh=mysh, pi=0.995, geno=geno, map=map, pheno=mypheno,        
            chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)  
}


### prediction
source("~/Documents/Rcodes/sh4Predict.R")
map="/home/NSF-SAM-GS/GenSel/rc_runs/SAM_map.txt"
for(i in 1:10){
  mysh <- paste("SAMrc_predict_s", i, ".sh", sep="")
  mypheno <- paste("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_s", i, "_val.txt", sep="")
  mymrkRes <- paste("/home/NSF-SAM-GS/GenSel/rc_runs/SAMrc_run41000_s", i, ".mrkRes1", sep="")
  sh4Predict(pwd=pwd, sh=mysh, geno=geno, pheno=mypheno, map=map, mrkRes=mymrkRes) 
                        
}


