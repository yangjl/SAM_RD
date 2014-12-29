### Jinliang Yang
### transform dsf to GenSel format!

test <- read.table("~/Documents/SAM_GS/snpmatrix/merged_chrall.dsf", header=TRUE, nrow=5)

#### recoding to GenSel format
python ~/Documents/PyCodes/dsf2GWAS/dsf2GWAS_v1.0.py -h
python ~/Documents/PyCodes/dsf2GWAS/dsf2GWAS_v1.0.py -i merged_chrall.dsf -o merged_chrall_samid.gensel -s 6 -e 374 -m 0


