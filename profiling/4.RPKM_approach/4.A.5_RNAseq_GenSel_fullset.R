# Jinliang Yang
# Purpose: GenSel of cross-validation 80% training and 20% val
# date: July.14.2014
# location: server.9


pwd="/home/NSF-SAM-GS/GenSel/rc_runs"
geno1="/home/NSF-SAM-GS/GenSel/rc_runs/SAM_trc_rpkm.txt.newbin"
geno2="/home/NSF-SAM-GS/GenSel/rc_runs/SAM_trc_rpkm_rescaled.txt.newbin"
map="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_map.txt"
mypheno <- "/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt"

source("~/Documents/Rcodes/sh4GenSel.R")
### BayesC run using raw geno
mysh <- "SAMrc41000_sam382_raw.sh"
sh4GenSel(pwd=pwd, sh=mysh, pi=0.995, geno=geno1, map=map, 
          pheno=mypheno, findsale ="no",       
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)  

### BayesC run using rescaled geno
mysh <- "SAMrc41000_sam382_rescaled.sh"
sh4GenSel(pwd=pwd, sh=mysh, pi=0.995, geno=geno2, map=map, 
          pheno=mypheno, findsale ="no",       
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)  


### add addMapInfoToMarkers yes
source("~/Documents/Rcodes/sh4GenSel.R")
mysh <- "SAMrc41000_sam382_rescaled_addmap.sh"
map="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_map.txt"
sh4GenSel(pwd=pwd, sh=mysh, pi=0.99, geno=geno2, map=map, 
          pheno=mypheno, findsale ="no",       
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)  





