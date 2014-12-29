# Jinliang Yang
# August 12nd, 2014

source("~/Documents/Rcodes/sh4GenSel.R")
### test runs:
sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/testrun", 
          sh="SAM_test.sh",
          geno="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.gensel", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.map",
          chainLength=1000, burnin=100, varGenotypic=0.2, varResidual=0.2)

### test runs: pi=0.995
sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/run_mergedSNPs/", 
          sh="SAM_run41000_pi995.sh", pi=0.995,
          geno="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.1)

sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/run_mergedSNPs/", 
          sh="SAM_run41000_pi9995.sh", pi=0.9995,
          geno="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.1)

sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/run_mergedSNPs/", 
          sh="SAM_run41000_pi999.sh", pi=0.999,
          geno="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.1)

#### addmap 10Mb
sh4GenSel(pwd="/home/NSF-SAM-GS/GenSel/run_mergedSNPs/", 
          sh="SAM_run41000_pi9995_win10m.sh", pi=0.9995,
          geno="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.gensel.newbin", 
          pheno="/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt",
          map="/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid_10m.map",
          chainLength=41000, burnin=1000, varGenotypic=0.2, varResidual=0.1)
