library(ggplot2)
library(data.table)
library(ggsci)
#source("./basics_TP.R")
source("Data/rec1EE/rec1EE_utils.R")


load("./Data/rec1EE/genotype/Rhom_poolseq.RData")



###################################################
#### EEV1403 haplotype around rec-1 

eev1403snps = fread(paste0("./Data/rec1EE/genotype/", "SMR_inbredlines_snps_CeMEEv2_ws245.csv"))
eev1403gt=fread(paste0("./Data/rec1EE/genotype/", "SMR_inbredlines_genotype_CeMEEv2_ws245.csv"))[,"EEV1401"]

#load("./Data/rec1_cemee_lines/riails_parental/rec1_riails_parental_genotype.Rdata")
mm = match(paste0(snps$chrom,snps$pos),paste0(eev1403snps$chrom, eev1403snps$pos))
eev1403gt = eev1403gt[mm,]
eev1403gt = as.numeric(unlist(c(eev1403gt)))

w.mut = grepl("LR", dimnames(ref)[[2]])
w.wt = grepl("HR", dimnames(ref)[[2]])
freqG0.wt = apply(alt[,w.wt,"G0"], 1, sum, na.rm=T) / (apply(ref[,w.wt,"G0"], 1, sum, na.rm=T) + apply(alt[,w.wt,"G0"], 1, sum, na.rm=T))
freqG0.mut = apply(alt[,w.mut,"G0"], 1, sum, na.rm=T) / (apply(ref[,w.mut,"G0"], 1, sum, na.rm=T) + apply(alt[,w.mut,"G0"], 1, sum, na.rm=T))
#freqG0.wt = apply(alt[,w.wt,"g0"], 1, sum, na.rm=T) / (apply(ref[,w.wt,"g0"], 1, sum, na.rm=T) + apply(alt[,w.wt,"g0"], 1, sum, na.rm=T))
#freqG0.mut = apply(alt[,w.mut,"g0"], 1, sum, na.rm=T) / (apply(ref[,w.mut,"g0"], 1, sum, na.rm=T) + apply(alt[,w.mut,"g0"], 1, sum, na.rm=T))

freqG0.mut.polarized = freqG0.mut 
freqG0.mut.polarized[which(eev1403gt==0)] = 1 - freqG0.mut.polarized[which(eev1403gt==0)]
freqG0.mut.polarized[is.na(eev1403gt) | eev1403gt==0.5] = NA

freqG0.wt.polarized = freqG0.wt 
freqG0.wt.polarized[which(eev1403gt==0)] = 1 - freqG0.wt.polarized[which(eev1403gt==0)]
freqG0.wt.polarized[is.na(eev1403gt)| eev1403gt==0.5] = NA


div = data.frame(pos = snps$pos, chrom = snps$chrom, freq.wt = freqG0.wt.polarized, freq.mut = freqG0.mut.polarized)
div$freq.all = (div$freq.wt + div$freq.mut)/2

#ggplot(div, aes(freq.all))+geom_histogram()
#ggplot(div, aes(pos, freq.all))+geom_point(size=0.1, alpha=0.1)+facet_wrap(~chrom)

#The minimum frequency of eev1403 allele we would expect from 
#sampling 340lines with an expected proportion of eev1403 haplotype of 0.5
#extreme = qbinom(0.001, size=340,prob = 0.5)/340

#error = which(div$freq.all<extreme) # likely genotyping error of eev1403
#div$freq.mut[error]=NA
#div$freq.wt[error]=NA
#div$freq.all[error]=NA
div$diff = as.numeric(div$freq.mut) - as.numeric(div$freq.wt)

div$diff


div$Drel = unlist(lapply(1:nrow(div), function(i){ 
  
  
  p12 = 0.5*div$freq.mut[i]
  p1 = 0.5
  p2 = div$freq.all[i]
  q1 = 1-p1
  q2 = 1-p2
  
  D = p12 - p1*p2
  Dmax = NA
  if(!is.na(D)){
    if(D>0){Dmax = min(p1*q2,q1*p2)}
    if(D<0){Dmax = min(p1*p2, q1*q2)}
  }
  
  D/Dmax
  
  #x = try(cor(x, wtx, use = "complete.obs"), silent = T)
  #if(class(x)=="try-error") x = NA
  #x
  
}))



##################################################################
####### PLOT DIFFERENCE IN MUTANT AND WT ANCESTRAL ALLELE FREQ ##
##################################################################

#ggplot(div, aes(pos/1e6, diff))+
#  geom_point( shape=1, size=0.2, alpha=0.5)+geom_hline(yintercept = 0, color="blue", linetype='dashed')+
#  geom_vline(data = data.frame(chrom='I', recpos=recpos), aes(xintercept = recpos/1e6), color='red')+
#  ylim(-1,1)+facet_grid(~chrom, scale="free")+
#  theme_minimal()+
#  xlab("Physical position (Mb)")+ylab("Difference in ancestral allele frequencies")


#ggplot(div, aes(diff))+
#  geom_histogram(bins=25)+facet_grid(~chrom, scale="free")+
#  theme_minimal()+
#  xlab("Physical position (Mb)")+ylab("Difference in ancestral allele frequencies")



################################################
### PLOT DIFFERENCE IN ALLELE FREQUENCIES


## MAIN PLOT: AF diff ~ pos


#ggplot(div, aes(pos/1e6, diff))+facet_grid(~chrom, scale="free_x", space="free_x")+
#  geom_point( shape=1, size=0.5/.pt, alpha=0.5)+geom_hline(yintercept = 0, color="blue", linetype='dashed')+
#  geom_vline(data = data.frame(chrom='I', recpos=recpos), aes(xintercept = recpos/1e6), color='red')+
#  theme_Publication3()+
#  facet_grid(~chrom, scale="free_x", space="free_x")+
#  scale_y_continuous(expand = c(0,0), lim=c(-1.2,1))+
#  xlab("Physical position (Mb)")+ylab("Allele frequency difference")






#pdiversity = ggplot(subset(div, freq.all<0.90), aes(pos/1e6, A6140.proportion))+
#  geom_point( shape=1, size=0.3)+geom_hline(yintercept = 0.5, color="blue", linetype='dashed')+
#  geom_vline(data = data.frame(chrom='I', recpos=recpos), aes(xintercept = recpos/1e6), color='red')+
#  ylim(0,1)+facet_grid(~chrom, scale="free_x", space = 'free_x')+
#  theme_Publication2()+
#  xlab("Physical position (Mb)")+ylab("A6140 haplotypes proportion")




#pD = ggplot(subset(div, freq.all<0.90), aes(pos/1e6, Drel))+
#  geom_point( shape=1, size=0.5/.pt, alpha=0.5)+
#  geom_hline(yintercept = 0, color="blue", linetype='dashed')+
#  geom_vline(data = data.frame(chrom='I', recpos=recpos), aes(xintercept = recpos/1e6), color='red')+
#  facet_grid(~chrom, scale="free_x", space = 'free_x')+
#  theme_Publication3()+
#  xlab("Physical position (Mb)")+ylab("D/Dmax")+ylim(-1,1)




#p=grid.arrange(pdiff+xlab(""),
#             #pdiversity+xlab("")+theme(strip.text.x = element_blank()),
#             pD+theme(strip.text.x = element_blank()), nrow=2, heights=c(1, 0.85))

#ggsave(filename="/Users/tomparee/Desktop/ancestral_diversity.png", plot=p, height=3.5, width=7.5, dpi=300)
#ggsave(filename="/Users/tomparee/Desktop/ancestral_diversity.pdf", plot=p, height=3.5, width=7.5, dpi=300)




mdiv=do.call(rbind,lapply(split(div, div$chrom), function(x){
  
  #x=split(div, div$chrom)[[1]]
  ppos = seq(0,max(x$pos), 5e3)
  #ppos1 = ppos-5e5
  #ppos2 = ppos+5e5
  
  meandiff=sapply(ppos, function(p){mean(subset(x, pos> p-5e5 & pos<=p+5e5)$diff, na.rm=T)})
  
  
  data.frame(chrom=x$chrom[1],
             pos = ppos,
             diff=meandiff)
  
}))



pDiff=ggplot(div, aes(pos/1e6, diff))+
  geom_point(data=data.frame(x=0,y=0), aes(x,y), color=NA)+
  geom_point( shape=1, size=0.3/.pt, alpha=0.3, color="grey40")+
  geom_hline(yintercept = 0, size = 1/.pt, linetype = 'dashed')+
  geom_point(data=data.frame(chrom="I", pos = recpos/1e6), aes(pos, y=1), shape=25, size=3/.pt, fill='black')+
  geom_text(data=data.frame(chrom="I", pos = recpos/1e6), aes(pos, y=1.2, label = "rec-1"),fontface="italic", size = 7/.pt, color='black')+
  theme_Publication3()+
  facet_grid(~chrom, scale="free_x", space="free_x")+
  scale_x_continuous(expand = c(0, 0),breaks = seq(0,25,4))+
  scale_y_continuous(expand = c(0, 0),breaks = c(-1,-0.5,0,0.5,1))+
  coord_cartesian(clip="off", ylim = c(-1.6,1))+
  xlab("Physical position (Mb)")+ylab("\u{0394} Allele frequency")+
  geom_line(data=mdiv, aes(pos/1e6,diff, color=chrom))+
  scale_color_npg()+theme_Publication3()+theme(legend.position="none")


## PLOT WITHIN MAIN PLOT: DENSITY OF AF diff
ds = list()

for(chr in 1:6){
  
  if(chr == 1){ br = c(-0.3, 0, 0.3, 0.6) }else{ br = c(-0.2, 0, 0.2) }
  d = ggplot(subset(div, chrom %in% numtorom(chr)), aes(diff))+
    geom_density(fill='lightgrey', size=1/.pt)+
    theme_Publication3()+
    geom_vline(xintercept = 0, color='black', linetype='dashed')+
    xlab("\u{0394} Allele freq.")+
    ylab("Density")+
    scale_x_continuous(breaks = br, expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,20,10))+
    theme(axis.text=element_text(size=4.5),
          axis.title=element_text(size=5,face="bold"),
          axis.text.x = element_text(angle = 90))
  
  
  
  
  ds[[chr]] = d
}



########################################
## ADD DENSITY PLOT WITHIN MAIN PLOT ###

library(grid)

p =pDiff + 
  annotation_custom2(grob=ggplotGrob(ds[[1]]), xmin = 2, xmax = 14.75,
                     data = data.frame(chrom='I', diff=1, pos=1),
                     ymin = -1.6, ymax = -0.35)+
  annotation_custom2(grob=ggplotGrob(ds[[2]]), xmin = 2, xmax = 14.75,
                     data = data.frame(chrom='II', diff=1, pos=1),
                     ymin = -1.6, ymax = -0.35)+
  annotation_custom2(grob=ggplotGrob(ds[[3]]), xmin = 0, xmax = 12.75,
                     data = data.frame(chrom='III', diff=1, pos=1),
                     ymin = -1.6, ymax = -0.35)+
  annotation_custom2(grob=ggplotGrob(ds[[4]]), xmin = 5.5, xmax = 17.25,
                     data = data.frame(chrom='IV', diff=1, pos=1),
                     ymin = -1.6, ymax = -0.35)+
  annotation_custom2(grob=ggplotGrob(ds[[5]]), xmin = 7, xmax = 19.75,
                     data = data.frame(chrom='V', diff=1, pos=1),
                     ymin = -1.6, ymax = -0.35)+
  annotation_custom2(grob=ggplotGrob(ds[[6]]), xmin = 4.5, xmax = 17.25,
                     data = data.frame(chrom='X', diff=1, pos=1),
                     ymin = -1.6, ymax = -0.35)







ggsave(filename="/Users/tomparee/Desktop/Figures/AllelefreqDiff.png", plot=p, height=1.9, width=5.4, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/Figures/AllelefreqDiff.pdf", plot=p,height=1.9, width=5.4, dpi=600)




mDrel=do.call(rbind,lapply(split(div, div$chrom), function(x){
  
  print(x$chrom[1])
  ppos = seq(0,max(x$pos), 5e3)
  #ppos1 = ppos-5e5
  #ppos2 = ppos+5e5
  
  meanD=sapply(ppos, function(p){mean(subset(x, pos> p-5e5 & pos<=p+5e5 & freq.all<0.90)$Drel, na.rm=T)})
  
  
  data.frame(chrom=x$chrom[1],
             pos = ppos,
             Drel=meanD)
  
}))

pD = ggplot(subset(div, freq.all<0.90), aes(pos/1e6, Drel))+
  geom_point(data=data.frame(x=0,y=0), aes(x,y), color=NA)+
  geom_point( shape=1, size=0.3/.pt, alpha=0.3, color="grey40")+
  geom_hline(yintercept = 0, size = 1/.pt, linetype = 'dashed')+
  geom_point(data=data.frame(chrom="I", pos = recpos/1e6), aes(pos, y=1.1), shape=25, size=3/.pt, fill='black')+
  geom_text(data=data.frame(chrom="I", pos = recpos/1e6), aes(pos, y=1.3, label = "rec-1"),fontface=4, size = 7/.pt, color='black')+
  theme_Publication3()+
  facet_grid(~chrom, scale="free_x", space="free_x")+
  scale_x_continuous(expand = c(0, 0),breaks = seq(0,25,4))+
  scale_y_continuous(expand = c(0, 0),breaks = c(-1,-0.5,0,0.5,1))+
  coord_cartesian(clip="off", ylim = c(-1,1))+
  xlab("Physical position (Mb)")+ylab("D'")+
  scale_color_npg()+theme(legend.position="none")+
  geom_line(data=mDrel, aes(pos/1e6,Drel, color=chrom))



ggsave(filename="/Users/tomparee/Desktop/ancestralD.png", plot=pD, height=1.22, width=5.4, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/ancestralD.pdf", plot=pD,height=1.22, width=5.4, dpi=600)





# 
# ################################################################################################################
# #### simulation tentative to estimate the expected difference of allele freq between HR and LR populations #####
# ################################################################################################################
# 
# # Warnings: it only estimate the difference in the F2 of the cross design
# # After that, lines were expanded from each single F2, so big opportunity for drift
# # + expanded lines were merged in equal proportions, but there is an error on worm number estimation
# # then populations were maintained a few generations to increases male frequencies, even though the population size were large
# # => So these simulations largely underestimate the expected allele frequency difference between HR and LR,
# #    as they lack a lot of steps. 
# 
# 
# 
# cemeesnps <- fread('./cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_RIL_snps_ws245.csv.gz')
# cemeegt <- fread('./cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_RIL_geno.csv.gz')
# cemeegt = as.matrix(cemeegt)
# 
# #Keep only lines from A6140 pop from which HR and lR populations were derived
# cemeegt= cemeegt[,grepl('A6140',colnames(cemeegt))]
# 
# eev1403snps = fread(paste0("./Data/rec1EE/genotype/", "SMR_inbredlines_snps_CeMEEv2_ws245.csv"))
# eev1403gt=fread(paste0("./Data/rec1EE/genotype/", "SMR_inbredlines_genotype_CeMEEv2_ws245.csv"))[,"EEV1401"]
# 
# 
# # #Keep only common snps
# mm = match(paste0(cemeesnps$chrom,cemeesnps$pos),paste0(eev1403snps$chrom, eev1403snps$pos))
# eev1403gt = eev1403gt[mm,]
# 
# cemeegt = cbind(eev1403gt,cemeegt)
# cemeegt = as.data.frame(cemeegt)
# cemeegtPol = cemeegt
# 
# cemeegtPol[which(cemeegtPol[,1]==0),] = 1 - cemeegtPol[which(cemeegtPol[,1]==0),]
# 
# ncross = 340
# chr = 'I'
# 
# 
# cemeegtchr = cemeegtPol[cemeesnps$chrom==chr,]
# cemeesnpschr = cemeesnps[cemeesnps$chrom==chr,]
# 
# F1 = data.frame(F1.g1 = rep(1,ncross),
#            F1.g2 = sample(2:ncol(cemeegt),ncross, replace=T))
# 
# 
# 
# simu = do.call(cbind, lapply(1:nrow(F1), function(i){
#   
#   f1 = F1[i,]
#   
#   geno = do.call(cbind, lapply(rep(c("wt","mut"), 2), function(recgeno){
#     
# 
#     recombinant = sample(c(T,F), 1)
#     
#     if(recgeno=='mut' & chr == 'I'){
#       whichgenome = 1
#     }else{
#       whichgenome = sample(1:2,1)
#     }
#     
#     
#     if(recombinant==T){
#       breakpoint = which.min(abs(cemeesnpschr$cM - runif(1)*50))
#       
#       g = c(cemeegtchr[1:breakpoint, unlist(c(f1[whichgenome]))],
#         cemeegtchr[(breakpoint+1):nrow(cemeegtchr), unlist(c(f1[1:2 != whichgenome]))])
#       
#     }else{
#       
#       g = cemeegtchr[,unlist(c(f1[whichgenome]))]
#     }
#     
#     
#   }))
#   
#   colnames(geno) = paste0(c("wt.g1_","mut.g1_","wt.g2_","mut.g2_"), i)
#   geno
# 
# }))
# 
# AFdiff = apply(simu[,grepl("mut",colnames(simu))], 1, mean,na.rm=T)-apply(simu[,grepl("wt",colnames(simu))], 1, mean,na.rm=T)
# hist(AFdiff)



#ggplot(cbind(cemeesnpschr, AFdiff), aes(pos, AFdiff))+geom_point()
# 
# ### Import snps ID and genotypes of the cemee
# cemeesnps <- fread('./cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_RIL_snps_ws245.csv.gz')
# cemeegt <- fread('./cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_RIL_geno.csv.gz')
# cemeegt = as.matrix(cemeegt)
# 
# #Keep only lines from A6140 pop from which HR and lR populations were derived
# cemeegt= cemeegt[,grepl('A6140',colnames(cemeegt))]
# 
# #Keep only common snps
# mm = match(paste0(cemeesnps$chrom,cemeesnps$pos),paste0(snps$chrom, snps$pos))
# 
# ###########################################################################
# #### % of polymorphism in HR and LR populations recovered from the cemee
# 1 - (sum(is.na(mm)) / length(mm))
# # 0.908553
# 
# 
# ###########################################################################
# #### Correlations between AFC in A6140 and HR/LR pops
# 
# refG0 = ref[,,"G0"]
# altG0 = alt[,,"G0"]
# freqG0 = apply(altG0, 1, sum, na.rm=T) / (apply(refG0, 1, sum, na.rm=T) + apply(altG0, 1, sum, na.rm=T))
# 
# cemeefreq = apply(cemeegt, 1, sum, na.rm=T) / apply(cemeegt, 1, function(x){ sum(!is.na(x))})
# 
# mm = match(paste0(cemeesnps$chrom,cemeesnps$pos),paste0(snps$chrom, snps$pos))
# freqG0 = freqG0[mm]
# 
# freqG0 = cbind(freqG0, cemeefreq)
# cor.test(freqG0[,1], freqG0[,2], use = "complete.obs")
# 
# 
# ###########################################################################
# #### Compare Major allele in A6140 and allele of EEV1403
# 
# load("./Data/rec1_cemee_lines/riails_parental/rec1_riails_parental_genotype.Rdata")
# mm = match(paste0(cemeesnps$chrom,cemeesnps$pos),paste0(snps$chrom, snps$pos))
# eev1403gt = geno[mm,]
# eev1403gt[is.na(mm),]=0 #would be better to use freebayes to genotype at these pos
# eev1403gt = eev1403gt[,"EEV1401"]
# eev1403gt = as.numeric(eev1403gt)
# 
# sum(round(freqG0[,"cemeefreq"]) == eev1403gt, na.rm=T) / sum(!is.na(round(freqG0[,"cemeefreq"])))
# 
# 
# 
# 
# ###########################################################################
# #### Compare HR and LR populations
# 
# w.mut = grepl("LR", dimnames(ref)[[2]])
# w.wt = grepl("HR", dimnames(ref)[[2]])
# 
# freq.wt = apply(alt[,w.wt,], 1, sum, na.rm=T) / (apply(ref[,w.wt,], 1, sum, na.rm=T) + apply(alt[,w.wt,], 1, sum, na.rm=T))
# freq.mut = apply(alt[,w.mut,], 1, sum, na.rm=T) / (apply(ref[,w.mut,], 1, sum, na.rm=T) + apply(alt[,w.mut,], 1, sum, na.rm=T))
# 
# ### Number of polymorphism not present in both populations 
# sum((freq.wt ==0 &  freq.mut != 0) | (freq.wt !=0 &  freq.mut == 0))
# 
# ### Difference at G0 ###
# freqG0.wt = apply(alt[,w.wt,"G0"], 1, sum, na.rm=T) / (apply(ref[,w.wt,"G0"], 1, sum, na.rm=T) + apply(alt[,w.wt,"G0"], 1, sum, na.rm=T))
# freqG0.mut = apply(alt[,w.mut,"G0"], 1, sum, na.rm=T) / (apply(ref[,w.mut,"G0"], 1, sum, na.rm=T) + apply(alt[,w.mut,"G0"], 1, sum, na.rm=T))
# 
# # Correlation
# cor.test(freqG0.wt, freqG0.mut)
# 
# # Average difference in AF
# mean(abs(freqG0.wt-freqG0.mut))
# sd(abs(freqG0.wt-freqG0.mut))
# 
