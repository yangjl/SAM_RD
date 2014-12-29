### Jinliang Yang
### 8/12/2014
### Xianran's merged SNP sets

### input genotype of SAM from Xianran merged dataset
setwd("~/Documents/SAM_GS/snpmatrix")
tab5rows <- read.table("merged_ch1.hmp", nrow=5, header=TRUE)
classes <- sapply(tab5rows, class)

#####
id_lookup <- function(infile=tab5rows){
  ## updated the name_lookup table
  nmtable <- read.csv("~/Documents/SAM_GS/SAM_proj/data/SAM_name_lookup.csv")
  
  gidtable <- data.frame(gid=names(infile), order=1:ncol(infile))
  nmtable <- nmtable[!is.na(nmtable$accession),]
  gidtable2 <- merge(gidtable, nmtable[, c("SAMID", "accession")], 
                     by.x="gid", by.y="accession", all.x=TRUE)
  
  gidtable2[is.na(gidtable2$SAMID), ]$SAMID <- gidtable2[is.na(gidtable2$SAMID), ]$gid
  gidtable2 <- gidtable2[order(gidtable2$order),]
  return(gidtable2)
}

gid <- id_lookup(infile=tab5rows)


###############
hmp2dsf <- function(infile="merged_ch1.hmp", gid=gid, outfile="merged_chrall.dsf", write=FALSE, ...){
  ### input genotype of SAM from Xianran merged dataset
  #setwd("~/Documents/SAM_GS/snpmatrix")
  tab5rows <- read.table(infile, nrow=5, header=TRUE)
  
  
  if(sum(names(tab5rows) %in% gid$gid) != nrow(gid)){
    stop("names not match!")
  }
  classes <- sapply(tab5rows, class)
  hmp <- read.table(infile, header=TRUE, colClasses=classes)
  names(hmp) <- gid$SAMID
  names(hmp)[1:8] <- c("snpid", "alleles", "chr", "pos", "major", "minor", "MAF", "missing")
  
  hmp$snpid <- paste(hmp$chr, hmp$pos, sep="_")
  hmp$major <- gsub("\\/.*", "", hmp$alleles)
  hmp$minor <- gsub("^.*\\/", "", hmp$alleles)
  hmp$MAF <- hmp$panelLSID
  hmp$missing <- hmp$QCcode
  
  if(write){
    write.table(hmp[, c(1,5:8,12:380)], outfile, row.names=FALSE, quote=FALSE, sep="\t", ...)
    message(sprintf("[ %s ] SNPs write to [ %s ]", nrow(hmp), outfile))
  }
  return(hmp[,1:9])
}

##############
setwd("~/Documents/SAM_GS/snpmatrix")
hmp1 <- hmp2dsf(infile="merged_ch1.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE)
hmp2 <- hmp2dsf(infile="merged_ch2.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp3 <- hmp2dsf(infile="merged_ch3.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp4 <- hmp2dsf(infile="merged_ch4.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp5 <- hmp2dsf(infile="merged_ch5.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp6 <- hmp2dsf(infile="merged_ch6.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp7 <- hmp2dsf(infile="merged_ch7.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp8 <- hmp2dsf(infile="merged_ch8.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp9 <- hmp2dsf(infile="merged_ch9.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)
hmp10 <- hmp2dsf(infile="merged_ch10.hmp", gid=gid, outfile="merged_chrall.dsf", write=TRUE,
                col.names=FALSE, append=TRUE)

map <- rbind(hmp1, hmp2, hmp3, hmp4, hmp5, hmp6, hmp7, hmp8, hmp9, hmp10)
names(map)[c(3,4,9)] <- c("AGPv2_chr", "AGPv2_pos", "type") 
write.table(map, "merged_chrall_samid.map", sep="\t", row.names=FALSE, quote=FALSE)

####
map <- read.table("/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid.map", header=T)
map$AGPv2_pos <- round(map$AGPv2_pos/10, 0)
write.table(map, "/mnt/02/yangjl/Documents/SAM_GS/snpmatrix/merged_chrall_samid_10m.map", 
            sep="\t", row.names=FALSE, quote=FALSE)


