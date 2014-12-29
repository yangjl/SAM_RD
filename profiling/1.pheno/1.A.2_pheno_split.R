### Jinliang Yang
### updated: 7/14/2014
### split into training and validation sets


splitSets <- function(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s1"){
  pn <- round(nrow(sam) * ptraining, 0)
  idx <- sample(1:nrow(sam), pn, replace=FALSE)
  train <- sam[idx, ]
  val <- sam[-idx,]
  
  ### output file:
  out1 <- paste(outfile, "_train.txt", sep="")
  out2 <- paste(outfile, "_val.txt", sep="")
  write.table(train, out1, sep="\t", row.names=FALSE, quote=FALSE)
  write.table(val, out2, sep="\t", row.names=FALSE, quote=FALSE)
}

### output phenotype files for GenSel
sam <- read.table("~/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt", header=TRUE)
set.seed(1234)
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s1")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s2")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s3")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s4")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s5")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s6")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s7")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s8")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s9")
splitSets(sam=sam, ptraining=0.8, outfile="~/Documents/SAM_GS/SAM_proj/reports/SAM_s10")


