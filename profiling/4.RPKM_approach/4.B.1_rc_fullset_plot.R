# Jinliang Yang
# Purpose: GenSel of cross-validation 80% training and 20% val
# date: July.24.2014
# location: server.9

getBayes <- function(pwd="/home/NSF-SAM-GS/GenSel/rc_runs",
                     inputfile="SAMrc41000_sam382_raw.mrkRes1",
                     map="SAM_map.txt"){
  setwd(pwd)
  res <- read.table(inputfile, header=TRUE)
  map <- read.table(map, header=TRUE)
  res2 <- merge(map, res[c("Genotype", "Effect", "ModelFreq")], by="Genotype")
  message(sprintf("input [ %s ]", nrow(res2)))
  
  names(res2)[2:3] <- c("chr", "pos")
  res2$chr <- as.numeric(as.character(sub("chr", "", res2$chr)))
  res2 <- subset(res2, !is.na(chr))
  message(sprintf("remove unknow chr, remainning [ %s ]", nrow(res2)))
  return(res2)
}

#############
source("~/Documents/Rcodes/quickMHTplot.R")
par(mfrow=c(1,1))
#### pai=0.995
bayes1 <- getBayes(pwd="/home/NSF-SAM-GS/GenSel/rc_runs",
                   inputfile="SAMrc41000_sam382_raw.mrkRes1",
                   map="SAM_map.txt")
#bayes1s <- subset(bayes1, ModelFreq > 0.005)
quickMHTplot(res=bayes1, main="with raw count", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")

#### pai=0.995
bayes2 <- getBayes(pwd="/home/NSF-SAM-GS/GenSel/rc_runs",
                   inputfile="SAMrc41000_sam382_rescaled.mrkRes1",
                   map="SAM_map.txt")
#bayes1s <- subset(bayes1, ModelFreq > 0.005)
quickMHTplot(res=bayes2, main="re-scale geno", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")
abline(h=0.08, lty=2, col="red")

bayes2 <- bayes2[order(bayes2$ModelFreq, decreasing=TRUE),]

######################
prepheno <- function(ghat="/home/NSF-SAM-GS/GenSel/SAM_run41000.ghatREL1",
                     cgrRes="/home/NSF-SAM-GS/GenSel/SAM_run41000.cgrRes1", ...){
  h <- read.table(ghat, skip=1)
  #PEV=var(g/y)
  names(h) <- c("Genotype", "gHat", "log2vol", "Fix", "meanBias", "PEV", "R2")
  p <- read.table(cgrRes, skip=2)
  h$gHat <- h$gHat + p$V2
  
  plot(h$gHat, h$log2vol, xlab="Predicted", ylab="Observed", ...)
  print(cor.test(h$gHat, h$log2vol))
  return(h)
}

### 

setwd("/home/NSF-SAM-GS/GenSel/rc_runs")
par(mfrow=c(1,2))
t1 <- prepheno(ghat="SAMrc41000_sam382_rescaled.ghatREL1",
               cgrRes="SAMrc41000_sam382_rescaled.cgrRes1", 
               main="rescaled geno")
t2 <- prepheno(ghat="SAMrc41000_sam382_raw.ghatREL1",
               cgrRes="SAMrc41000_sam382_raw.cgrRes1", 
               main="raw geno")

gene2 <- subset(bayes2, ModelFreq > 0.08)
gene2 <- gene2[order(gene2$chr, gene2$pos),]

merge(gene1, gene2, by.x="Gene", by.y="Genotype")

des <- read.delim("~/db/AGPv2/ZmB73_5a_gene_descriptors.txt", header=TRUE)
bayes2a <- merge(bayes2, des, by.x="Genotype", by.y="gene_name", all.x=TRUE)
bayes2a <- bayes2a[order(bayes2a$ModelFreq, decreasing=TRUE),]
write.table(bayes2a, "~/Documents/SAM_GS/SAM_proj/reports/samvolume_rcmethod.txt", sep="\t",
            row.names=FALSE, quote=FALSE)
