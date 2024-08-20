##############################################################################
##################### afc r2 ~ delta r #######################################
library(data.table)
library(ggplot2)

PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

nsnp_window = 1000 # window size in number of snp to calculate ld between pairwise SNP combinations in each genomic region
#(if too big, it will take a long time to run, pairwise => increases exponentially)
#(if too small, sample only very close SNP)
nboot = 1000

load("genotype/Rhom_poolseq.RData")


freq = alt/(ref+alt)



#########################################################################
############# FUNCTIONS #################################################
#########################################################################

calculate.AFcor = function(freq, nboot = 10){
  #freq = matrix[snp, samples]
  cormatrix = cor(t(freq) , use = "pairwise.complete.obs")
  cormatrix=cormatrix^2
  cormatrix[lower.tri(cormatrix, diag=T)]=NA
  
  data = data.frame(r2 = mean(cormatrix, na.rm=T), nobs = sum(!is.na(cormatrix)))
  
  bootstraped.mean = unlist(lapply(1:nboot, function(b){
    is = sort(sample(1:nrow(freq), nrow(freq), replace=T))
    bootcor = cormatrix[is,is]
    mean(mean(bootcor, na.rm=T))
    }))
  
  return(list(data, bootstraped.mean))
}



#########################################################################
#############    RUN & SAVE    ##########################################
#########################################################################


afr2 = do.call(rbind, lapply(c("HR", "LR"), function(poptype){
  
  do.call(rbind, lapply(unique(snps$chrom), function(chr){
    
    print(chr)
    
    output = do.call(rbind, lapply(c(-1,1,0), function(signDeltaR){
      
      print(signDeltaR)
      
      if(signDeltaR==1) ix = snps$rrdiff > 0 & snps$chrom == chr
      if(signDeltaR==-1) ix = snps$rrdiff < 0 & snps$chrom == chr
      if(signDeltaR==0){ix = snps$chrom == chr; signDeltaR = "genome-wide"}
      
      freqx = freq[ix,grepl(poptype, dimnames(freq)[[2]]),]
      freqx = arraytomatrix(freqx, c(3,2))
      freqx = freqx[, apply(freqx, 2, function(x){sum(is.na(x))}) < nrow(freqx)]
      windows = get.win(size=nrow(freqx), winsize=nsnp_window, minsize = nsnp_window/2)
      
      AFcor =  do.call(rbind, lapply(windows, function(win){
        #win = windows[[1]]
        freqw = freqx[win,]
        data = calculate.AFcor(freq=freqw, nboot = nboot)
        data = unlist(data)
      }))
      
     
      AFcor = as.data.frame(AFcor)
      
      AFcor = c(nobs = sum(AFcor[,2]) , apply(AFcor[,-2],2, function(x){weighted.mean(x, AFcor[,2])}))
      AFcor = as.data.frame(matrix(AFcor, nrow=1))
      colnames(AFcor) = c("nobs", "mean.r2", paste0("boot", 1:(ncol(AFcor)-2)) )
      AFcor = cbind(chrom = chr, rrchange = signDeltaR , rec=poptype, AFcor)
      return(AFcor)
    }))
    
    output
    
  }))
  
  
}))





afr2$rec = ifelse(afr2$rec=="HR", "wt", "mut")

#save(afr2, file = "observed_r2/afr2.Rdata")
#write.table(afr2, file = "observed_r2/afr2.txt")
#load(file = "observed_r2/afr2.Rdata")

