PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
library(data.table)
library(ggplot2)
library(gridExtra)
source("./rec1EE_utils.R")

load("recmaps.Rdata")

# Calculate recombination rates from recombination maps
recrates = do.call(rbind,lapply(split(recmaps, paste0(recmaps$rec, recmaps$chrom)), function(x){
  #x = split(recmaps, paste0(recmaps$cross, recmaps$chrom))[[1]]
  pdis = diff(x$pos)
  rr = diff(x$genetic)/pdis
  pos1=x$pos[1:(nrow(x)-1)]
  pos2=x$pos[2:nrow(x)]
  
  out = data.frame(pos1 = pos1, pos2=pos2, rr=rr)
  out$chrom = x$chrom[1]
  out$cross = x$cross[1]
  out$rec = x$rec[1]
  out$background = x$background[1]
  out
  
}))

recrates = recrates[order(recrates$rec, recrates$pos1),]

# Calculate delta r (difference recombination rates mutant - wild-type)
recchange = recrates[recrates$rec == 'wt',]
colnames(recchange)[colnames(recchange) == "rr"] = "rr_wt"
recchange$rr_mut = recrates$rr[recrates$rec == 'mut']
recchange$rrdiff = recchange$rr_mut-recchange$rr_wt
recchange$rrdiff=recchange$rrdiff*1e6 #cM per Mb

save(recchange, file = "deltaR.Rdata")




