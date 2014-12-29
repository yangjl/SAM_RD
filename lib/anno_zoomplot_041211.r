#! /usr/bin/Rscript --vanilla 
# --default-packages=utils,stats,lattice,grid,getopts
# need to check if the line above works on the web deployment machine.
# Modified by Jinliang Yang
# 4.10.2011

#require(stats);
#require(utils);
require(grid);
#require(lattice);
#install.packages(getopt);
#require(getopts);
# server: 129.186.85.42
setwd("~/GWAS/plot")


### DATA Formatting:
####################################################################################
###GWAS by PLINK
plink <- read.table("wmean_merged_GWAS.qassoc", header=T)
dim(plink)
#[1] 844065      9
plink$ID <- paste(plink$CHR, plink$BP, sep="_")

### modelFreq from GenSel
Freq <- read.table("KRN_wmean_HapApex_half.mrkRes1", header=T)
tail(Freq)
Freq$ID <- paste(Freq$chr, Freq$pos, sep="_")

pval <- merge(plink, Freq, by="ID")

### founder
founder <- read.delim("HapMap1_Apex2_040711_founders.txt", header=TRUE)

pval <- merge(pval, founder, by.x="SNP", by.y="rs")
dim(pval)
#[1] 844063     52
pval <- pval[,c(1:10,13:22,25:52)]
pval <- pval[order(pval$CHR, pval$BP),]
#set1: HapMap1 set2: Apex, set12: both
write.table(pval, "pval_master.csv", sep=",", row.names=FALSE,quote=FALSE)
####################################################################################
pval <- read.csv("pval_master.csv")
# get 30 random SNP
idx <- sample(1:844063, 30)
control <- pval[idx,]





####################################################################################
# Gene Model: B73_V1 (1st isoform from Sanzhen):
fgs <- read.table("~/GWAS/Anno/ZmB73_4a.53_filtered_genes.first.transcript.gff", header=TRUE)
dim(fgs)
#[1] 493914     12
exon <- fgs[fgs$Feature=="exon" | fgs$Feature=="intron",]

######################################################################
library(ggplot2)
exonlen <- fgs[fgs$Feature=="exon",]
exonlen$dist<- exonlen$End-exonlen$Start+1

res <- ddply(exonlen, .(Gene), summarise,
	len = sum(dist))

res$Gene <- gsub("_T01", "",res$Gene)
res$Gene <- gsub("_T02", "",res$Gene)
res$Gene <- gsub("_T03", "",res$Gene)
res$Gene <- gsub("_T04", "",res$Gene)
res$Gene <- gsub("_T05", "",res$Gene)
res$Gene <- gsub("_T06", "",res$Gene)
res$Gene <- gsub("_T07", "",res$Gene)
res$Gene <- gsub("_T08", "",res$Gene)
res$Gene <- gsub("_T11", "",res$Gene)
res$Gene <- gsub("_FGT", "_FG",res$Gene)

rna <- read.csv("/home/yangjl/GWAS/Anno/B73_ear+tassel-AGPv1_FGSv4a-2.csv")
dim(rna)
#[1] 32540

###################################################################
gene <- fgs[fgs$Feature=="gene",]
gene <- merge(gene, rna, by.x="Gene", by.y="FGS_v4a")
gene <- gene[,c(1:5,10:12)]
names(gene) <- c("Gene","Chr","Feature","Start","End", "Strand","Ear","Tassel")
gene <- gene[order(gene$Chr, gene$Start),]
genebkup <- gene
gene <- genebkup

gene <- merge(gene, res, by="Gene", all=TRUE)
dim(gene)

	sum1 <- sum(gene$Ear)
	sum2 <- sum(gene$Tassel)
	gene$Ear <- round((1000000000*gene$Ear)/(gene$len*sum1),0)
	gene$Tassel <- round((1000000000*gene$Tassel)/(gene$len*sum2),0)
gene <- gene[order(gene$Chr, gene$Start),]
write.table(gene, "FGS_V1_Ear_Tassel_RPKM.csv", sep=",", row.names=FALSE)
####################################################################################
gene <- read.csv("FGS_V1_Ear_Tassel_RPKM.csv")













#########################################################################
pdf(file="GWAS_Zoomplot_v1.pdf", paper="a4r")
########### Chr1
# QTL 1
zoomplot(mb =150000, chr=1, snppos=7210988)
zoomplot(mb =150000, chr=1, snppos=7532782)
#zoomplot(mb =150000, chr=1, snppos=23803429)
zoomplot(mb =100000, chr=1, snppos=12214500)
zoomplot(mb =150000, chr=1, snppos=23663730)


# QTL 2
zoomplot(mb =150000, chr=1, snppos=165957057)
zoomplot(mb =250000, chr=1, snppos=167192700)
zoomplot(mb =150000, chr=1, snppos=179743994)

# QTL 3
zoomplot(mb =150000, chr=1, snppos=223809294)

# QTL 4
zoomplot(mb =150000, chr=1, snppos=281550408)
zoomplot(mb =150000, chr=1, snppos=281789566)
zoomplot(mb =150000, chr=1, snppos=281931181)
zoomplot(mb =150000, chr=1, snppos=282263931)

# QTL 0
zoomplot(mb =150000, chr=1, snppos=201990240)
zoomplot(mb =150000, chr=1, snppos=202226724)

########### Chr2
# QTL 5
zoomplot(mb =50000, chr=2, snppos=15912082)
zoomplot(mb =150000, chr=2, snppos=16650439)
zoomplot(mb =150000, chr=2, snppos=22543227)

# QTL 0

zoomplot(mb =500000, chr=2, snppos=94981521)
zoomplot(mb =250000, chr=2, snppos=213175830)
zoomplot(mb =150000, chr=2, snppos=223399678)
zoomplot(mb =150000, chr=2, snppos=223605304)
zoomplot(mb =150000, chr=2, snppos=225334100)

########### Chr3
# QTL 6
zoomplot(mb =150000, chr=3, snppos=159310805)
zoomplot(mb =50000, chr=3, snppos=165439145, myylim=c(0,0.3))

# QTL 7
zoomplot(mb =50000, chr=3, snppos=192266158)

# QTL 8
zoomplot(mb =150000, chr=3, snppos=6885330, myylim=c(0,0.3))
zoomplot(mb =150000, chr=3, snppos=20108088, myylim=c(0,0.3))
zoomplot(mb =150000, chr=3, snppos=25509654, myylim=c(0,0.3))


########### Chr4
# QTL 9
zoomplot(mb =50000, chr=4, snppos=3397611)
zoomplot(mb =150000, chr=4, snppos=4223353)

# QTL 10
zoomplot(mb =100000, chr=4, snppos=5301191)
zoomplot(mb =150000, chr=4, snppos=40640711)
zoomplot(mb =150000, chr=4, snppos=41174823)
#zoomplot(mb =150000, chr=4, snppos=122529884)
zoomplot(mb =150000, chr=4, snppos=203305659)
zoomplot(mb =150000, chr=4, snppos=204074816)
zoomplot(mb =150000, chr=4, snppos=204687128, myylim=c(0,0.8))#most significant one
zoomplot(mb =150000, chr=4, snppos=204764297)
zoomplot(mb =150000, chr=4, snppos=207014128)
zoomplot(mb =150000, chr=4, snppos=227245183)
zoomplot(mb =150000, chr=4, snppos=241373338,myylim=c(0,0.8))
zoomplot(mb =150000, chr=4, snppos=242139950)

########### Chr5
zoomplot(mb =150000, chr=5, snppos=18480132)
zoomplot(mb =150000, chr=5, snppos=93378634)
zoomplot(mb =150000, chr=5, snppos=105502729)
zoomplot(mb =150000, chr=5, snppos=211860602)
zoomplot(mb =50000, chr=5, snppos=214249467)

########### Chr6
zoomplot(mb =150000, chr=6, snppos=80388440)
zoomplot(mb =150000, chr=6, snppos=122546064)
zoomplot(mb =150000, chr=6, snppos=128737704)
zoomplot(mb =100000, chr=6, snppos=150901962)

########### Chr7
zoomplot(mb =150000, chr=7, snppos=20599112)
zoomplot(mb =150000, chr=7, snppos=88043386)
zoomplot(mb =150000, chr=7, snppos=100192900)
zoomplot(mb =100000, chr=7, snppos=153896451)

########### Chr8
zoomplot(mb =150000, chr=8, snppos=117027428)
#zoomplot(mb =150000, chr=8, snppos=151830409)
zoomplot(mb =150000, chr=8, snppos=151830434)
#zoomplot(mb =100000, chr=8, snppos=164401324)
#zoomplot(mb =150000, chr=8, snppos=164408442)
zoomplot(mb =100000, chr=8, snppos=164408470)

########### Chr9
zoomplot(mb =150000, chr=9, snppos=3776348)
zoomplot(mb =150000, chr=9, snppos=4924167)
zoomplot(mb =150000, chr=9, snppos=96351970)
zoomplot(mb =150000, chr=9, snppos=148836165)
zoomplot(mb =100000, chr=9, snppos=150066752)

########### Chr10
zoomplot(mb =150000, chr=10, snppos=18720422)
zoomplot(mb =150000, chr=10, snppos=19388294)
zoomplot(mb =150000, chr=10, snppos=23784833)
zoomplot(mb =150000, chr=10, snppos=88599962)
zoomplot(mb =150000, chr=10, snppos=121185188)
zoomplot(mb =150000, chr=10, snppos=137386351)
zoomplot(mb =150000, chr=10, snppos=138519111)
zoomplot(mb =150000, chr=10, snppos=138951633)

dev.off()




























