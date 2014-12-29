# Jinliang Yang
# updated: 7/14/2014
# Purpose: 5-fold cross-validation

#65%
####################
get_predict <- function(ghat="/home/NSF-SAM-GS/GenSel/train_val/SAM_predict_s10.ghat1",
                        cgrRes="/home/NSF-SAM-GS/GenSel/train_val/SAM_run41000_s10.cgrRes1", ...){
  h <- read.table(ghat, skip=1)
  #PEV=var(g/y)
  names(h) <- c("Genotype", "gHat", "log2vol", "Fix")
  p <- read.table(cgrRes, skip=2)
  h$gHat <- h$gHat + p$V2
  
  r <- round(cor(h$gHat, h$log2vol), 2);
  
  plot(h$gHat, h$log2vol, main=paste("r =", r), xlab="Predicted", ylab="Observed", ...)
  # plot 
  abline(lm(h$log2vol ~ h$gHat), col="red", lwd=2)
  print(cor.test(h$gHat, h$log2vol))
  return(h)
}

### 
par(mfrow=c(3,3))
for(i in 2:10){
  myghat <- paste("/home/NSF-SAM-GS/GenSel/rc_runs/SAMrc_predict_s", i, ".ghat2", sep="")
  mycgr <- paste("/home/NSF-SAM-GS/GenSel/rc_runs/SAMrc_run41000_s", i, ".cgrRes1", sep="")
  tem <- get_predict(ghat=myghat, cgrRes=mycgr, type="p", pch=16, col="blue")
}
