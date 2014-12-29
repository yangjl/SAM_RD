# Jinliang Yang
# updated: 7/24/14
# AGPv2 Annotation

######################## AGPv2 Chr length ###################################
AGPv2_chrlength <- read.table("~/Documents/sharedata/ZmB73_RefGen_v2.length", header=FALSE)
AGPv2_chrlength <- AGPv2_chrlength
AGPv2_chrlength <- t(AGPv2_chrlength);

######################## Filtering the duplicated annotation #################
FGSv2 <- read.table("~/db/AGPv2/ZmB73_5b_FGS.gff", header=FALSE)
names(FGSv2) <- c("seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attribute")
FGSv2_exon <- subset(FGSv2, feature=="exon" & seqname!= "UNKNOWN")
write.table(FGSv2_exon, "~/db/AGPv2/ZmB73_5b_FGS_exon.gff", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

FGSv2_info <- read.table("~/db/AGPv2/ZmB73_5b_FGS_info.txt", header=TRUE)

######################## Working Gene Set #################
WGSv2 <- read.table("~/db/AGPv2/ZmB73_5a_WGS.gff", header=FALSE)
names(WGSv2) <- c("seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attribute")
WGSv2_exon <- subset(WGSv2, feature=="exon" & seqname!= "UNKNOWN")
write.table(WGSv2_exon, "~/db/AGPv2/ZmB73_5a_WGS_exon.gff", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

WGSv2_info <- read.delim("~/db/AGPv2/ZmB73_5a_WGS_info.txt", header=TRUE)

######################## 5a named genes #################
namedgenes <- read.delim("~/db/AGPv2/ZmB73_5a_named_genes.txt", header=TRUE)
dim(namedgenes)
#[1] 620   6

#source("~/Documents/Rcodes/save.append.R")
#save.append(list=c("AGPv2_chrlength"), file="test.RData",
#            description="first load")
