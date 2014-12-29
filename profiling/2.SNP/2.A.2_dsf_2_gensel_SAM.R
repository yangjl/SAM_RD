### Jinliang Yang
### transform dsf to GenSel format!

#### using genotype id to generate F1 ids
seudoF1 <- function(ids = saml$SAMID){
  
  f1 <- data.frame(p1=rep(ids, each=45), p2=rep(ids, times=45))
  f1$o1 <- as.numeric(as.character(gsub("SAM", "", f1$p1)))
  f1$o2 <- as.numeric(as.character(gsub("SAM", "", f1$p2)))
  f1 <- f1[order(f1$o1, f1$o2),]
  
  f1 <- subset(f1, o1 < o2)
  f1$genotype <- paste(f1$p1, f1$p2, sep="x")
  return(f1[, c(1,2,5)])
}


######## 
sam2 <- read.csv("~/Documents/SAM_GS/SAM_proj/data/SAM_large_small.csv")
saml <- subset(sam2, Instruction == "large")
#A679 is duplicated
saml <- saml[!duplicated(saml$SAMID),] 
#[1] 45  3
genol <- seudoF1(ids=saml$SAMID)
write.table(genol, "~/Documents/SAM_GS/SAM_proj/reports/SeudoF1_large.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

##### impute density SNPs
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.0.py -d SeudoF1_large.txt -i \
SAM_SNPs_chrall_samid.dsf -o SAM_large_chrall.gensel -s 6 -e 385 --header yes -m 0



######### SAM small
sams <- subset(sam2, Instruction == "small")
#[1] 45  3
genos <- seudoF1(ids=sams$SAMID)
write.table(genos, "~/Documents/SAM_GS/SAM_proj/reports/SeudoF1_small.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

##### impute density SNPs
python ~/Documents/PyCodes/impute4diallel/impute4diallel_v1.0.py -d SeudoF1_small.txt -i \
SAM_SNPs_chrall_samid.dsf -o SAM_small_chrall.gensel -s 6 -e 385 --header yes -m 0


