### Jinliang Yang
### August 27th, 2014

nms <- c("win", "start", "end", "snps", "var", "cumvar", 
         "p0", "pavg", "map_pos1", "map_pos2", "chr_mb")
winrc <- read.table("/home/NSF-SAM-GS/GenSel/rc_runs/SAMrc41000_sam382_rescaled_addmap.winQTL2",
                    header=FALSE, skip=1)
names(winrc) <- nms
### heritability 0.60

winsnp <- read.table("/home/NSF-SAM-GS/GenSel/run_mergedSNPs/SAM_run41000_pi9995_win10m.winQTL1",
                     header=FALSE, skip=1)
names(winsnp) <- nms
###heritability 0.74

wins <- merge(winrc[, c("var", "chr_mb")], winsnp[, c("var", "chr_mb")], by="chr_mb", all=TRUE)

plot(wins$var.x, wins$var.y, xlim=c(0, 1.8), ylim=c(0, 1.8), )
abline(0,1,col="red")
cor.test(wins$var.x, wins$var.y)



