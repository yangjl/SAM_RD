# Jinliang Yang
# July 10th, 2014

source("~/Documents/Rcodes/sh4GenSel.R")

traits <- c("arc_length", "area", "Height", "HtRa", "para_coeff", 
            "Radius", "surface_area", "volume")
varRes <- c(1623, 926207, 361, 0.016, 0.006, 36.5, 101700000, 1.50964e+11)
varGeno <- c(1623, 926207, 361, 0.016, 0.006, 36.5, 101700000, 1.50964e+11) 
################ testrun for add genotype
for(i in 1:8){
  mysh <- paste(traits[i], "_run41000.sh", sep="")
  mygeno <- "/home/NSF-SAM-GS/snpmatrix/merged_chrall_samid_maf1miss6.gensel.newbin"
  mypheno <- paste("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/reports/SAM_pheno382_", traits[i],
                   ".txt", sep="")
  
  sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/realrun", 
            sh=mysh, pi=0.9995, geno=mygeno, pheno=mypheno,
            map="/home/NSF-SAM-GS/snpmatrix/merged_chrall_samid_maf1miss6.map",
            chainLength=41000, burnin=1000, varGenotypic=varGeno[i], varResidual=varRes[i])
  
}




