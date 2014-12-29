### Jinliang Yang
### transform dsf to GenSel format!

test <- read.table("~/Documents/SAM_GS/snpmatrix/merged_chrall_maf1miss6.dsf", header=TRUE, nrow=5)

#### recoding to GenSel format
python ~/Documents/PyCodes/dsf2GWAS/dsf2GWAS_v1.0.py -h
python ~/Documents/PyCodes/dsf2GWAS/dsf2GWAS_v1.0.py -i merged_chrall_maf1miss6.dsf \
-o merged_chrall_samid_maf1miss6.gensel -s 6 -e 374 -m 0


