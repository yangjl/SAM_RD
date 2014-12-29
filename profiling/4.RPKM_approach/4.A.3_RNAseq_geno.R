# Jinliang Yang
# updated: 7/14/2014
# Purpose: quick plot of GWAS results

ob <- load("~/Documents/SAM_GS/SAM_proj/cache/SAM_rc_rpkm.RData")

names(nrc)[1] <- "Genotype"
trc <- t(nrc[, c(1,7:ncol(nrc)) ])
#dim(trc)
idx <- which(is.na(trc[,1]))
write.table(trc[-idx,], "~/Documents/SAM_GS/SAM_proj/cache/SAM_trc_rpkm.csv", sep=",",
            col.names=FALSE, row.names=TRUE, quote=FALSE)

################
# map file:
map <- rc[,1:3]
names(map) <- c("Genotype", "AGPv2_chr", "AGPv2_pos")
map$AGPv2_chr <- gsub("chr", "", map$AGPv2_chr)
map <- map[order(map$AGPv2_chr, map$AGPv2_pos),]
map$AGPv2_pos <- as.numeric(as.character(map$AGPv2_pos))
#map$Genotype <- gsub("\\.", "_", map$Genotype)
map$bin <- paste(map$AGPv2_chr, round(map$AGPv2_pos/1000000, 0), sep="_")
map$AGPv2_pos <- round(map$AGPv2_pos/10,0)
write.table(map, "~/Documents/SAM_GS/SAM_proj/reports/SAM_map.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

######################################################
#rescale it to -10 to 10;
rescaleRPKM <- function(){
  trc <- read.csv("~/Documents/SAM_GS/SAM_proj/cache/SAM_trc_rpkm.csv", header=TRUE)
  # remove non SAM id with NA
  #trc <- trc[-94,]
  write.table(trc, "~/Documents/SAM_GS/SAM_proj/reports/SAM_trc_rpkm.txt", sep="\t",
              row.names=FALSE, quote=FALSE)
  
  ### rescale to -10 -> 10
  source("~/Documents/Rcodes/rescale.R")
  trc2 <- trc;
  for(i in 2:ncol(trc2)){
    xrange <- range(trc2[,i])
    if(xrange[1] != xrange[2]){
      trc2[,i] <- rescale(trc2[,i], c(-10, 10))
    } 
  }
  write.table(trc2, "~/Documents/SAM_GS/SAM_proj/reports/SAM_trc_rpkm_rescaled.txt", sep="\t",
              row.names=FALSE, quote=FALSE)
}

rescaleRPKM()



