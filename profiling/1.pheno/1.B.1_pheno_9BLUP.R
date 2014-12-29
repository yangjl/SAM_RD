# Jinliang Yang
# created: June 9, 2014
# phenotype of SAM volume

sam <- read.csv("~/Documents/SAM_GS/SAM_proj/data/SAM_BLUP.csv")

#pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.F1_pheno.pdf", width=10, height=5)
#layout(matrix(c(1,1,2,3,4,1,1,5,6,7), 2, 5, byrow = TRUE))
par(mfrow=c(3,3))
#plotKRN()
plot(density(sam$Height), col="black", lwd=4, bty="n", main="Height", xlab="")
plot(density(sam$Radius), col="black", lwd=4, bty="n", main="Radius", xlab="")
plot(density(sam$arc_length), col="black", lwd=4, bty="n", main="arc_length", xlab="")
plot(density(sam$HtRa), col="black", lwd=4, bty="n", main="HtRa", xlab="")
plot(density(sam$diameter), col="black", lwd=4, bty="n", main="diameter", xlab="")
plot(density(sam$area), col="black", lwd=4, bty="n", main="area", xlab="")
plot(density(sam$volume), col="black", lwd=4, bty="n", main="volume", xlab="")
plot(density(sam$para_coeff), col="black", lwd=4, bty="n", main="para_coeff", xlab="")
plot(density(sam$surface_area), col="black", lwd=4, bty="n", main="surface_area", xlab="")


##### correlation plot
source("~/Documents/Rcodes/Correlation_plot.R")
# read in the data
#
#pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.F2_cor.pdf", height=8, width=8)
pairs(sam[, 2:10], text.panel = diag, upper.panel=panel.smooth, 
      lower.panel=panel.cor, gap=0, main="", pch=19, col="grey", lwd=2)


########

write_pheno <- function(sam1=sam[,-5], pwd="~"){
  nm <- read.csv("~/Documents/SAM_GS/SAM_proj/data/SAM_name_lookup.csv", header=T)
  sam2 <- merge(nm[,c("accession", "SAMID")], sam1, by.x="accession", by.y="ID")
  sam2 <- sam2[order(sam2$SAMID), ]
  for(i in 3:ncol(sam2)){
    ### output phenotype files for GenSel
    myout <- paste(pwd, "SAM_pheno382_", names(sam2)[i], ".txt", sep="")
    sam3 <- sam2[, c(2, i)]
    sam3$Fix <- 1
    names(sam3)[1] <- "Genotype"
    write.table(sam3, myout, sep="\t",
                row.names=FALSE, quote=FALSE)
  }
}

write_pheno(sam1=sam[,-5], pwd="~/Documents/SAM_GS/SAM_proj/reports/")
  
  
  



