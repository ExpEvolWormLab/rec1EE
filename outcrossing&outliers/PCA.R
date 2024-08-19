library(ggfortify)
library(data.table)
library(dplyr)
library(data.table)

path = "./Data/rec1EE/"

source(paste0(path,"utils_poolseq_analysis.R"))
load(paste0(path,"Rhom_poolseq.RData"))

freq = alt/(alt+ref)
coverage = alt + ref

rec1linked = snps$pos < 5e6 & snps$chrom == 'I'
lowcov = apply(coverage, 1, function(x) sum(x<10, na.rm=T))
lowcov = lowcov > 0

freq = freq[!rec1linked & !lowcov,,]



#logit = F
#if(logit){
#extreme = apply(freq, 1, function(x) sum(x==0 | x==1, na.rm=T))
#extreme = extreme > 0
#f0 = qlogis(freq[!extreme,,"g0"])
#f2 = qlogis(freq[!extreme,,g2])
#}


pca = do.call(rbind, lapply(c('g17', 'g24', 'g33', 'g41'), function(g2){
  
  print(g2)
  
  f0 = freq[,,"g0"] 
  f2 = freq[,,g2]
  
  # If sequencing is missing intrapolate
  miss = apply(freq[,,g2], 2, function(x) sum(!is.na(x))) == 0
  wg2 = which(dimnames(freq)[[3]] == g2)
  
  if(g2 != 'g41'){
    
    dgen = diff(as.numeric(tstrsplit(dimnames(freq)[[3]], "g")[[2]])[c(wg2-1, wg2, wg2+1)])
    p = dgen/sum(dgen)
    
    f2[,miss]=freq[,miss,wg2-1] * (1-p[1]) +  freq[,miss,wg2+1] * (1-p[2]) 
    
  }else{
    f2[,miss] = freq[,miss,'g33']
  }
  
  
  # Allele frequency change
  afc = t(f2-f0)
  afc = as.data.frame(afc)
  afc = cbind(assign.infopop(data.frame(pop=rownames(afc))), afc)
  
  # Take out SNVs unsequenced in at least one pop
  wna2 = apply(afc, 2 , function(x) sum(is.na(x)))
  wna2 = wna2 > 0
  afc = afc[,!wna2]
  
  df = do.call(rbind, lapply(c(T,F), function(cc){
    
    # Perform PCA: centered and uncentered
    afc.pca <- prcomp(afc[,4:(length(afc))], center = cc,scale = F)
    
    #summary(afc.pca)
    df.pca = data.frame(PC1 = afc.pca$x[,1]/(afc.pca$sdev[1] * sqrt(nrow(afc))),
                        PC2 = afc.pca$x[,2]/(afc.pca$sdev[2] * sqrt(nrow(afc))))
    df.pca$pop = afc$pop
    df.pca$rec_env = paste0(afc$rec, afc$env)
    
    if(cc){df.pca$center = 'centered'}else{df.pca$center = 'uncentered'}
    df.pca 
    
  }))
  
  df$generation = g2
  
  df
}))






pca$rec_env = factor(pca$rec_env, labels = c('Mutant x Domestication',
                                             'Mutant x Novel',
                                             'Wild-type x Domestication',
                                             'Wild-type x Novel'))



source("./basics_TP.R")
p=ggplot()+
  theme_Publication3()+
  geom_point(data = subset(pca, center == "centered"), aes(PC1, PC2, color = rec_env), size=3/.pt)+
  geom_text(data = subset(pca, pop %in% c('LR1', 'SHR4') & center == "centered"), 
            aes(PC1, PC2+0.2, color = rec_env, label = pop),size=7/.pt,show.legend = F)+
  scale_color_manual(values=c('skyblue', blueWT, "#F4BB6C",orangeMut), 
                     breaks=c('Wild-type x Domestication', 'Wild-type x Novel', 'Mutant x Domestication', 'Mutant x Novel'), name="")+
  xlim(c(min(pca$PC1)-0.1, max(pca$PC1)+0.1))+
  facet_wrap(~toupper(generation), scales = 'free')+
  theme(panel.border=element_rect(linetype=1, fill=NA))+
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        strip.text = element_text(size = 7))+
  ylim(range(subset(pca, center == "centered")$PC2)*1.3)+
  xlim(range(subset(pca, center == "centered")$PC1)*1.3)


ggsave(filename="/Users/tomparee/Desktop/PCA.png", plot=p, width=3.5,height=2.2, dpi=500)
ggsave(filename="/Users/tomparee/Desktop/PCA.pdf", plot=p,  width=3.5,height=2.2, dpi=500)

  



