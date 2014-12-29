### Jinliang Yang
### July 21, 2014

####################################################################################
ZoomInPlot <- function(pval=pval, gene=gene, chr="chr4", window=500000, 
                       snppos=165957057, ylim=c(0,0.5), idxmax=NULL){
  ## pval: pval of the file, CHR, BP
  ## gene: gene model
  
  startpos = snppos - window
  endpos = snppos + window
  startmb = startpos/1000000
  endmb = endpos/1000000
  
  pval <- pval[pval$chr == chr & pval$pos >= startpos & pval$pos <= endpos,]
  mygene <- gene[gene$chr == chr & gene$start >= startpos & gene$end <= endpos,]
  
  ######
  grid.newpage();
  ####plot layout
  mylayout=grid.layout(5,3, 
                       widths=unit(c(4,1,2),c('lines','null','line')),
                       heights=unit(c(3, 3, 1, 1, 5), 
                                    c('lines', 'null','lines','null','lines'))
  )
  pushViewport(viewport(layout=mylayout));
  #### 
  ## ModelFreq----------------------------------------
  BayesVp=dataViewport(
    xscale=c(startmb, endmb), yscale=ylim,
    extension=c(0,.05),
    layout.pos.row=2,layout.pos.col=2,
    name="bayes",
    clip="off");
  
  pushViewport(BayesVp);
  #grid.xaxis();
  grid.yaxis();
  grid.rect();
  grid.points(x=pval$pos/1000000,y=pval$ModelFreq, pch=19, 
              gp=gpar(col="red", alpha=0.5, cex=0.5));	
  
  grid.lines(x=unit(c(0,1), "npc"), y=0.02,
             gp=gpar(col="grey", lwd=2, lty=2),
             default.units='native');
  #### highlight the most significant one!
  if(!is.null(idxmax)){
    grid.points(x=pval$pos[idxmax]/1000000,y=pval$ModelFreq[idxmax], gp=gpar(col="red"));	
    grid.text(label=pval$snpid[idxmax], x=pval$pos[idxmax]/1000000, 
              y=pval$ModelFreq[idxmax]+0.1, default.units="native")
  }
  upViewport()
  
  #### y-axis
  vp41=viewport(layout.pos.row=2,layout.pos.col=1, name="vp41",clip="off");
  pushViewport(vp41);
  grid.text(label="Model Frequency",x=unit(1, 'lines'), just="center", rot=90)
  upViewport()
  
  ## Gene Model
  GeneVp=dataViewport(xscale=c(startmb, endmb), yscale=c(0, 2), extension=c(0,.05),
                      layout.pos.row=4, layout.pos.col=2, name="gene", clip="off");
  pushViewport(GeneVp);
  grid.xaxis();
  #grid.yaxis();
  grid.rect();
  
  for(i in 1:nrow(mygene)){
    grid.lines(x=c(mygene$start[i]/1000000, mygene$end[i]/1000000), y=(i%%2) + 0.4,
               gp=gpar(col="darkblue", fill= "darkblue", lwd=3, lineend="square"),
               default.units='native');
    grid.text(as.character(mygene$Gene[i]), x=(mygene$start[i]+mygene$end[i])/2000000, 
              y=(i%%2) + 0.7, default.units="native",
              check.overlap=TRUE, gp=gpar(cex=0.7));
  }
  upViewport()
  
  pushViewport(viewport(layout.pos.row=5, layout.pos.col=2, clip="off"));
  grid.text(label=paste("Position on ", chr, " (Mb)", sep=""), just="center")
  upViewport(0)
  
  return(mygene)
}# End of the function


