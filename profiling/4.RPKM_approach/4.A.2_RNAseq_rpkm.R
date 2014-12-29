# Jinliang Yang
# updated: 7/14/2014
# Purpose: quick plot of GWAS results

###################################
norRPKM <- function(rc=rc){
  exonsize <- read.csv("~/Documents/XP-GWAS/XP-GWAS_proj/cache/FGSv2_canonical_exonSize.csv")
  rc <- merge(exonsize[, c("gene_id", "exonSize")], rc, by.x="gene_id", by.y="Gene")
  source("~/Documents/XP-GWAS/XP-GWAS_proj/lib/get_rpkm.R")
  
  rc <- get_rpkm(data=rc, colrange=7:ncol(rc))
  message(sprintf("[ %s ] genes were RPKM normalized!", nrow(rc)))
  return(rc)
}

############
rc <- read.table("~/Documents/SAM_GS/snpmatrix/NSF-SAM.380.individuals.FGSv5b.60.readcounts.txt",
                 header=TRUE)

nrc <- norRPKM(rc=rc)

###################################
id_lookup <- function(nrc=nrc){
  ## updated the name_lookup table
  nmtable <- read.csv("~/Documents/SAM_GS/SAM_proj/data/SAM_name_lookup.csv")
  
  gidtable <- data.frame(gid=names(nrc), order=1:ncol(nrc))
  nmtable <- nmtable[!is.na(nmtable$ID_in_snpfile),]
  gidtable2 <- merge(gidtable, nmtable[, c("SAMID", "ID_in_snpfile")], 
                     by.x="gid", by.y="ID_in_snpfile", all.x=TRUE)
  
  gidtable2[is.na(gidtable2$SAMID), ]$SAMID <- gidtable2[is.na(gidtable2$SAMID), ]$gid
  gidtable2 <- gidtable2[order(gidtable2$order),]
  return(gidtable2)
}

gid <- id_lookup(nrc=nrc)
names(nrc) <- gid$SAMID

##########
source("~/Documents/Rcodes/save.append.R")
save.append(list=c("rc", "nrc"), file="~/Documents/SAM_GS/SAM_proj/cache/SAM_rc_rpkm.RData",
            description=c("original read count", "normalized read count with SAM ID!"))