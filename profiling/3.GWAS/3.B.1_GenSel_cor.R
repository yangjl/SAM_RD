## Jinliang Yang
## 7/24/2014
## correlation of training and validation

#65%

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

setwd("/home/NSF-SAM-GS/GenSel")
par(mfrow=c(1,3))
t1 <- prepheno(ghat="SAM_run41000.ghatREL1",
               cgrRes="SAM_run41000.cgrRes1", 
               main="pai=0.995 (3,500 SNPs)")

t3 <- prepheno(ghat="SAM_run41000_pi.ghatREL2",
               cgrRes="SAM_run41000_pi.cgrRes2", 
               main="pai=0.999 (700 SNPs)")

t2 <- prepheno(ghat="SAM_run41000_pi.ghatREL1",
               cgrRes="SAM_run41000_pi.cgrRes1", 
               main="pai=0.9995 (350 SNPs)")


