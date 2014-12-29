getBayes <- function(inputfile="/home/NSF-SAM-GS/GenSel/SAM_run41000.mrkRes1"){
  res <- read.table(inputfile, header=TRUE)
  message(sprintf("input [ %s ]", nrow(res)))
  res$chr <- gsub("_.*", "", res$snpid)
  res$chr <- as.numeric(as.character(sub("chr", "", res$chr)))
  res <- subset(res, !is.na(chr))
  res$pos <- as.numeric(as.character(gsub(".*_", "", res$snpid)))
  message(sprintf("remove unknow chr, remainning [ %s ]", nrow(res)))
  return(res)
}