### Jinliang Yang
### July 25th, 2014
### purpose: association of phenotype and read count!

RC1qtl <- function(pgdat, ...){
  ########### model part! #############
  # Y = tissue + RPKM
  # Note: col2=response, col3=explainary variable
  fit <- glm(pgdat[,2] ~ pgdat[,3], data=pgdat, ...);
  ##########################################
  #summary(fit)
  pval <- try(summary(fit)$coefficients[2,"Pr(>|t|)"], silent=TRUE);
  if(is.null(pval)){
    return(NA)
  }else{
    return(pval)
  }
}

###########
pheno <- read.table("~/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt", header=TRUE)
trc <- read.csv("~/Documents/SAM_GS/SAM_proj/cache/SAM_trc_rpkm.csv", header=TRUE)
ob <- load("~/Documents/SAM_GS/SAM_proj/cache/SAM_rc_rpkm.RData")
map <- nrc[, 1:6]

map$pval <- sapply(1:nrow(map), 
                   function(i){
                     pg <- merge(pheno[,1:2], trc[, c("Genotype", map$gene_id[i])], by="Genotype")
                     res <- RC1qtl(pg, family = gaussian)
                     return(res);
                   })

map$pval <- as.numeric(as.character(map$pval))
map$qval <- p.adjust(map$pval, method="bonferroni")
hist(map$pval)
map <- map[order(map$pval),]

nrow(subset(map, qval<=0.01))
mapbkup <- map

map$log10p <- -log10(map$pval)
names(map)[3:4] <- c("chr", "pos")
map <- map[!is.na(map$log10p),]
map$chr <- as.numeric(as.character(gsub("chr", "", map$chr)))
map <- map[!is.na(map$chr),]
#############
source("~/Documents/Rcodes/quickMHTplot.R")
par(mfrow=c(1,1))
#bayes1s <- subset(bayes1, ModelFreq > 0.005)
quickMHTplot(res=map, main="glm of read count", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="-log10p", yaxis=NULL,
             col2plot="log10p")
