# according to Dr. Nettleton's code
############################
get_chistat <- function(snpdf=snp50[, c("high_ref", "high_alt")]){
  #out <- snp50[, 1:5]
  ## high vs. low
  snpdf <- as.matrix(snpdf)
  stat <- apply(snpdf, 1, function(y){
    #Example data vector:
    #y=c(42,20,30,35,18,50)
    #Arranged in a table like those on the board last week (for illustration)
    #tab=matrix(y,byrow=T,ncol=2)
    #tab
    
    x=1:3
    tab=matrix(y,byrow=T,ncol=2)
    o=glm(cbind(tab[,1],tab[,2])~x,
          family=binomial(link=logit))
    chisqStat=anova(o)[2,2]
    chisqStat
  })  
  return(stat)
}
