# Jinliang Yang
# July 10th, 2014

source("~/Documents/Rcodes/sh4GenSel.R")
### test runs:
sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/testrun", 
          sh="SAM_test.sh",
          geno="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.gensel", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map",
          chainLength=1000, burnin=100, varGenotypic=1, varResidual=2)

### test runs: pi=0.995
sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/", 
          sh="SAM_run41000.sh", pi=0.995,
          geno="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)

sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/", 
          sh="SAM_run41000_pi.9995.sh", pi=0.9995,
          geno="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)

sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/", 
          sh="SAM_run41000_pi.999.sh", pi=0.999,
          geno="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)

#### addmap 10Mb
sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/", 
          sh="SAM_run41000_pi9995_win10m.sh", pi=0.9995,
          geno="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/home/NSF-SAM-GS/snpmatrix/SAM_SNPs_chrall_samid_10m.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.2)
