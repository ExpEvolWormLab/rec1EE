PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")



load("indirect_selection/MR_rec1freq.Rdata")

scoef = do.call(rbind, lapply(split(rec1freq, rec1freq$replicate), function(x){
  s = try(summary( glm(rec1freq~generation, data=x, family="binomial") )$coef[2,1], silent = T)
  afc = try(summary(lm(rec1freq~generation, data=x))$coef[2,1], silent = T)
  if(class(s)[1]=='try-error'){s = NA}
  if(class(afc)[1]=='try-error'){afc = NA}
  x
  x = x[1, -which(colnames(x) %in% c("generation", "rec1freq", "sample"))]
  x$scoef = as.numeric(s)
  x$afc =  as.numeric(afc)
  x
  
}))

save(scoef, file="indirect_selection/MR_selectionCoeficients.Rdata")


