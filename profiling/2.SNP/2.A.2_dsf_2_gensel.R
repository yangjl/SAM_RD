### Jinliang Yang
### transform dsf to GenSel format!

setwd("/home/NSF-SAM-GS/snpmatrix/")
test <- read.table("SAM_SNPs_chrall_samid.dsf", header=TRUE, nrow=5)

#### recoding to GenSel format
python ~/Documents/PyCodes/dsf2GWAS/dsf2GWAS_v1.0.py -h
python ~/Documents/PyCodes/dsf2GWAS/dsf2GWAS_v1.0.py -i SAM_SNPs_chrall_samid.dsf -o SAM_SNPs_chrall_samid.gensel -s 6 -e 385 -m 0
ll

