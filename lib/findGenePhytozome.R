#ap <- paste()
#Files in the annotation subdirectory:
#  1) Zmays_284_6a.annotation_info.txt
#A summary of annotation details available in Phytozome. This is a tab-delimited file, as follows:
#  (Note: Columns are blank if no corresponding data is available)
#1: Phytozome internal transcript ID (potentially useful to connect to biomart datasets)
#2: Phytozome gene locus name
#3: Phytozome transcript name
#4: Phytozome protein name (often same as transcript name, but this can vary)
#5: PFAM
#6: Panther
#7: KOG
#8: KEGG ec
#9: KEGG Orthology
#10: Gene Ontology terms (NOTE: these are automated results from interpro2go in most genomes, *not* empirically derived)
#11: best arabidopsis TAIR10 hit name
#12: best arabidopsis TAIR10 hit symbol
#13: best arabidopsis TAIR10 hit defline
#14: best rice hit name
#15: best rice hit symbol
#16: best rice hit defline

##############################################################################
findGenePhytozome <- function(fgs=NULL, sig=sig, binsize=100000, annotation=TRUE){
  ### Note chr should be [chr1, chr10...]
  if(is.null(fgs)){
    fgs <- read.table("~/bin/NGSbin/fgsv2.txt", header=TRUE);
    message("#FGSv2 loaded automatically!")
  }
  
  message(sprintf("#[ %s ] SNPs > 0.02 cutoff", nrow(sig)))
  out <- data.frame()
  for(i in 1:nrow(sig)){
    mys <- sig$pos[i]-binsize
    mye <- sig$pos[i]+binsize
    myg <- subset(fgs, chr==sig$chr[i] & start > mys & end < mye)
    if(nrow(myg) > 0 ){
      myg$SNP <- sig$snpid[i]
    }
    out <- rbind(out, myg)
  }
  out <- out[!duplicated(out$Gene),]
  message(sprintf("#[ %s ] unique genes were identified with binsize [ %s ]", nrow(out), binsize))
  out <- out[order(out$chr, out$start),]
  out$order <- 1:nrow(out)
  if(annotation==TRUE){
    des <- read.delim("~/db/Zmays_284_6a/6a/annotation/Zmays_284_6a.annotation_info.txt", header=F)
    names(des) <- c("PhytoID", "geneid", "txid", "proid", "PFAM", "Panther", "KOG", "KEGGec",
                    "KEGGorth", "GO", "ara_hit_id", "ara_hit_symbol", "ara_hit_defline",
                    "rice_hit_id", "rice_hit_symbol", "rice_hit_defile")
    
    
    #des$geneid <- gsub("_.*", "", des$gene_name)
    #xref <- read.delim("~/db/AGPv2/ZmB73_5a_xref.txt", header=TRUE)
    #xref$geneid <- gsub("_.*", "", xref$transcript)
    ###
    #test <- merge(des, xref, by="geneid")
    
    out <- merge(out, des, by.x="Gene", by.y="geneid", all.x=TRUE) 
  }
  message(sprintf("#[ %s ] Phytozome annotations [Zmays_284_6a] were found!", nrow(out)))
  return(out[order(out$order),])
  
}

