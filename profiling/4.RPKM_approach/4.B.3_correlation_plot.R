### Jinliang Yang
### 8/20/2014

##############
rawrpkm <- read.table("~/Documents/SAM_GS/SAM_proj/reports/SAM_trc_rpkm.txt", header=TRUE)
mypheno <- read.table("/mnt/02/yangjl/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt", header=TRUE)


genes ?



pdf("~/Documents/SAM_GS/SAM_proj/graphs/correls.pdf", width=8, height=8)
for(i in 1:nrow(gene2)){
  gp <- merge(mypheno, rawrpkm[, c("Genotype", gene2$Genotype[i])], by="Genotype")
  tit <- paste(gene2$Genotype[i], "\n", 
               "r =", round(cor(gp[,4], gp[,2]),2),"  ", 
               "p-value =", sprintf("%.1e", cor.test(gp[,4], gp[,2])$p.value)
  )
  x=gp[,4]
  y=gp[,2]
  plot(x, y, xlab="RPKM", ylab="log2(Volume of SAM)", main=tit)
  
  fit <- lm( y ~ x)
  abline(fit, col='red', lwd=2)
}
dev.off()

