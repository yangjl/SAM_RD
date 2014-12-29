# Jinliang Yang
# created: June 9, 2014
# phenotype of SAM volume

sam <- read.csv("~/Documents/SAM_GS/SAM_proj/data/SAM_methylsal_volume.csv")

par(mfrow=c(1,2))
hist(sam$mean_volume, breaks=30, main="mean volume", xlab="mean", col="lightgrey")
hist(log2(sam$mean_volume), breaks=30, main="log2 mean volume", xlab="log2(mean)", col="lightgrey")
# 382 5

#### normallity test 1: Shapiro's test,
shapiro.test(sam$mean_volume)
shapiro.test(log2(sam$mean_volume))

#install.packages("nortest")
library(nortest)
#### normality test 2: Anderson Darling
#http://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r
ad.test(sam$mean_volume)
ad.test(log2(sam$mean_volume))
qqnorm(log2(sam$mean_volume))
hist(log2(sam$mean_volume), breaks=30)

sam$log2vol <- log2(sam$mean_volume)

### output a name lookup table:
sam$order <- paste("SAM", sam$order, sep="")
write.table(sam, "~/Documents/SAM_GS/SAM_proj/data/SAM_name_lookup.csv", sep=",",
            row.names=FALSE, quote=FALSE)

### output phenotype files for GenSel
sam2 <- sam[, c("order", "log2vol")]
sam2$Fix <- 1
names(sam2)[1] <- "Genotype"
write.table(sam2, "~/Documents/SAM_GS/SAM_proj/data/SAM_pheno382.txt", sep="\t",
            row.names=FALSE, quote=FALSE)
