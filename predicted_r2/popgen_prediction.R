PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")
load("genotype/Rhom_poolseq.RData")
Rsnps = snps 

load("./recombinationMaps/recmaps.Rdata")


nLoci=5000
Ne <- 1e3
s <- 0.01 # only relevant for Roze 2021 equation
h <- 0.5  # only relevant for Roze 2021 equation

dbin = 0.001 # Linkage disequilibrium bin size
#dbins = rev(-seq(0,0.05,dbin))
dbins = seq(0,0.05,dbin)


### Caluclate Otha and Kimura sigma_d^2 = 1/4*Ne*c
### & Roze, 2021 estimation of negative LD
### # Hill robertson (1968 ? )  see doi: 10.1534/genetics.118.300642

interference = lapply(1:1000, function(b){
  
  if(b %% 10 == 0) print(b)
  OUT=try(do.call(rbind,lapply(c('I','II','III','IV', "V", "X"), function(chr){
    
    #print(chr)
    
    # Select n evenly spaced snps
    tarsnps = subset(Rsnps, chrom==chr)
    n=nLoci
    tarsnps = tarsnps[order(tarsnps$pos),]
    tarsnps = tarsnps[sort(sample(1:nrow(tarsnps), n)),] 
    
    do.call(rbind, lapply(c('wt', 'mut'), function(recgeno){
      
      #Calculate genetic distance from recombination maps
      map =subset(smoothmaps, chrom==romtonum(chr) & rec==recgeno)
      map = map[order(map$pos),]
      map$genetic = map$genetic*50
      
      gpos = approx(map$pos, map$genetic, xout = tarsnps$pos, rule=2)$y
      gdis = c(dist(gpos))
      r = gdis/100
      minr = 1e-2 
      r[r<minr]=minr
      
      D_otha <- mean(1/(4*Ne*r),na.rm=T) # Otha and Kimura
      Er2_HillRobertson <- mean(1/(1+(4*Ne*r)), na.rm=T) # Hill robertson (1968)  see doi: 10.1534/genetics.118.300642
      signed_D_Roze <- mean((s^2*h*(1-4*h))/(2*Ne*(r+2*s*h)^2*(r+3*s*h)),na.rm=T) # roze
      
      
      data.frame(chrom=chr, rec=recgeno,
                 od2_otha=D_otha,
                 Er2_HillRobertson=Er2_HillRobertson,
                 signed_D_Roze=signed_D_Roze,
                 nsampledSNP=n, Ne=Ne, s=s, h=h)
      
    }))
  })))
  
  if(class(OUT)=='try-error') return(NULL)
  
  OUT$nboot = b
  
  return(OUT)
})



interference=do.call(rbind, interference)



save(interference,file=paste0( "./Data/rec1EE/mixing_efficiency/predictedInterference_Ne",Ne, "_s",s,"_h", h,"_nloci", nLoci, ".Rdata" ))


# histrf=do.call(rbind, lapply(c('wt', 'mut'), function(recgeno){
#   
#   rf=unlist(lapply(c('I','II','III','IV', "V", "X"), function(chr){
#     
#     print(chr)
#     
#     # Select n evenly spaced snps
#     tarsnps = subset(Rsnps, chrom==chr)
#     n=1000
#     tarsnps = tarsnps[order(tarsnps$pos),]
#     tarsnps = tarsnps[floor(seq(1, nrow(tarsnps), nrow(tarsnps)/n)),] # n equally spaced SNV in cumulative SNV dist
#     
#     #Calculate genetic distance from smoothed maps
#     map =subset(smoothmaps, chrom==romtonum(chr) & rec==recgeno)
#     map = map[order(map$pos),]
#     map$genetic = map$genetic*50
#     
#     gpos = approx(map$pos, map$genetic, xout = tarsnps$pos, rule=2)$y
#     gdis = c(dist(gpos))
#     r = gdis/100
#     r
#   }))
#   
#   
#   rbin = 0.01
#   rbins = seq(0,0.5,rbin)
#   nperbins = unlist(lapply(1:(length(rbins)-1), function(i){
#     #print(i)
#     x1 = rbins[i]
#     x2 = rbins[i+1]
#     sum(rf>=x1 & rf<x2)
#   }))
#   
#   out = setNames(as.data.frame(cbind(rbins[1:(length(rbins)-1)],rbins[2:length(rbins)], nperbins)), c("rf1","rf2","n"))
#   out$rec=recgeno
#   out
#   
# })) 
# 
# 
# 
# 
# 
# 


##########################################################################
################## PLOT ##################################################
##########################################################################

meanHR=NULL
for(nLoci in c(10,100,1000)){
  load(file=paste0( "./Data/rec1EE/mixing_efficiency/predictedInterference_Ne1000_s0.01_h0.5_nloci", nLoci, ".Rdata" ))
  meanHR=rbind(meanHR,aggregate(Er2_HillRobertson ~ nboot+rec+nsampledSNP,interference, mean))
}


load(file=paste0( "./Data/rec1EE/mixing_efficiency/predictedInterference_Ne1000_s0.01_h0.5_nloci", 1000, ".Rdata" ))
meanHR=aggregate(Er2_HillRobertson ~ nboot+rec+nsampledSNP,interference, mean)
meanHR$rec=factor(meanHR$rec, levels = c('wt', 'mut'), labels = c("Wild-type", "Mutant"))

p=ggplot(meanHR, aes(as.factor(nsampledSNP), Er2_HillRobertson, color=rec))+
  stat_summary(fun.y = mean, 
               geom = "point", size=2/.pt, position = position_dodge(width=0.17))+ 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)*1.96, 
               fun.ymax = function(x) mean(x) + sd(x)*1.96, 
               geom = "errorbar",
               width=0.3/.pt,
               size=1/.pt,
               position = position_dodge(width=0.17))+
  scale_color_manual(values=c(blueWT, orangeMut),name="")+
  theme_Publication3()+
  labs(y=expression(bold(E(r^{2}))))+
  xlab("Number of loci under selection")#+
#theme(legend.position = "none")#+
#scale_y_continuous(expand = c(0.3,0))


#meanOtha=aggregate(od2_otha ~ nboot+rec,interference, mean)
#meanHR=aggregate(Er2_HillRobertson ~ nboot+rec,interference, mean)
#meanRoze=aggregate(signed_D_Roze ~ nboot+rec,interference, mean)

#library(ggbeeswarm)

meanHR$rec=factor(meanHR$rec, levels = c('wt', 'mut'), labels = c("Wild-type", "Mutant"))

p=ggplot(meanHR, aes(as.factor(nsampledSNP), Er2_HillRobertson, color=rec))+
  stat_summary(fun.y = mean, 
               geom = "point", size=2/.pt, position = position_dodge(width=0.17))+ 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)*1.96, 
               fun.ymax = function(x) mean(x) + sd(x)*1.96, 
               geom = "errorbar",
               width=0.3/.pt,
               size=1/.pt,
               position = position_dodge(width=0.17))+
  scale_color_manual(values=c(blueWT, orangeMut),name="")+
  theme_Publication3()+
  labs(y=expression(bold(E(r^{2}))))+
  xlab("Number of loci under selection")#+
  #theme(legend.position = "none")#+
  #scale_y_continuous(expand = c(0.3,0))


# ggplot(meanHR)+
#   geom_density(aes(Er2_HillRobertson, fill=rec), alpha=0.5, color='black')+
#   scale_color_manual(values=c(blueWT, orangeMut))+
#   theme_Publication3()+
#   labs(x=expression(bold(E(r^{2}))))+
#   xlab("")+theme(legend.position = "none", axis.text.x = element_text(face = "bold"))+
#   coord_cartesian(expand=0)+
#   scale_fill_manual(values=c(blueWT, orangeMut))
# 
# ggplot(meanHR, aes(factor(rec, levels = c('wt', 'mut')), Er2_HillRobertson, color=rec))+
#   stat_summary(fun.y = mean, 
#                geom = "point")+ 
#   stat_summary(fun.y = mean,
#                fun.ymin = function(x) mean(x) - sd(x)*1.96, 
#                fun.ymax = function(x) mean(x) + sd(x)*1.96, 
#                geom = "errorbar",
#                width=0.3)+
#   scale_color_manual(breaks = c("wt", "mut"),values=c(blueWT, orangeMut))+
#   theme_Publication3()+
#   labs(y=expression(bold(sigma[d]^{2})))+
#   xlab("")+theme(legend.position = "none")#+
# #scale_y_continuous(expand = c(0.3,0))


ggsave(filename="/Users/tomparee/Desktop/Figures/experctedr2.pdf", plot=p, height=1.3, width=2.7, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/Figures/experctedr2.png", plot=p, height=1.3, width=2.7, dpi=600)


ggplot(meanRoze, aes(factor(rec, levels = c('wt', 'mut')), signed_D_Roze*-1, color=rec))+
  stat_summary(fun.y = mean, 
               geom = "point")+ 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)*1.96, 
               fun.ymax = function(x) mean(x) + sd(x)*1.96, 
               geom = "errorbar",
               width=0.3)+
  scale_color_manual(breaks = c("wt", "mut"),values=c(blueWT, orangeMut))+
  theme_Publication3()+
  ylab("Negative LD")+
  xlab("")+theme(legend.position = "none")





prf=ggplot()+theme_Publication3()+
  geom_rect(aes(xmin=0,xmax=0.03, ymin=0, ymax=Inf), color="grey90", fill="grey90")+
  geom_step(data=subset(histrf, rec=="mut"),aes(x=rf1,y=n), alpha=1, color=orangeMut, size=0.5)+
  geom_rect(data=subset(histrf, rec=="mut"),aes(xmin=rf1,xmax=rf2, ymin=0, ymax=n),alpha=0.6, color=NA, fill=orangeMut)+
  geom_segment(data=subset(histrf, rec=="mut"), aes(x=rf1, xend=rf2, ,y=n, yend=n), alpha=1, color=orangeMut, size=0.5)+
  geom_step(data=subset(histrf, rec=="wt"),aes(x=rf1,y=n), alpha=1, color=blueWT, size=0.5)+
  geom_rect(data=subset(histrf, rec=="wt"),aes(xmin=rf1,xmax=rf2, ymin=0, ymax=n),alpha=0.6, color=NA, fill=blueWT)+
  geom_segment(data=subset(histrf, rec=="wt"), aes(x=rf1, xend=rf2, ,y=n, yend=n), alpha=1, color=blueWT, size=0.5)+
  coord_cartesian(expand=0)+
  scale_y_continuous(label=scientific_10)+
  xlab("Recombination probability")+ylab("Count")



nlowrecombing = aggregate(n~rec, subset(histrf, rf2<=0.03), sum)

plowrf = ggplot(nlowrecombing, aes(factor(rec, levels = c('wt', 'mut')), n))+
  theme_Publication3()+
  geom_bar(stat="identity", width=0.5, color="black", aes(fill=rec))+
  scale_fill_manual(breaks = c("wt", "mut"),values=c(blueWT, orangeMut))+
  theme(legend.position = "none")+
  scale_y_continuous(expand=c(0,0), label=scientific_10)+
  xlab("")+ylab("Count")
  



r <- seq(0.0001,0.5,0.0001)
negative_D_roze <- ((s^2*h*(1-4*h))/(2*Ne*(r+2*s*h)^2*(r+3*s*h)))*-1
D_otha <- 1/(4*Ne*r)
Drf = data.frame(rf=r, negative_D_roze=negative_D_roze, D_otha=D_otha)


proze = ggplot(subset(Drf, rf<=0.05), aes(rf, negative_D_roze ))+
  geom_rect(aes(xmin=0,xmax=0.03, ymin=0, ymax=Inf), color="grey90", fill="grey90")+
  geom_line(size=1)+theme_Publication3()+coord_cartesian(expand=0)+xlim(0,0.05)+
  labs(y="Negative LD", x="Recombination probability")
  

potha = ggplot(subset(Drf, rf<=0.05), aes(rf, D_otha ))+
  theme_Publication3()+
  geom_rect(aes(xmin=0,xmax=0.03, ymin=0, ymax=Inf), color="grey90", fill="grey90")+
  geom_line(size=1)+
  coord_cartesian(expand=0)+xlim(0,0.05)+
  scale_y_continuous(trans='sqrt')+
  labs(y=expression(bold(sigma[d]^{2})), x="Recombination probability")


ppopgen = gridExtra::grid.arrange(potha, proze, ncol=1)

p=gridExtra::grid.arrange(prf+theme(plot.margin = unit(c(5,20,5,5), "pt")),
                        plowrf+theme(plot.margin = unit(c(5,25,5,5), "pt")),
                        ppopgen, nrow=1, widths=c(1,0.4,0.5))




ggsave(filename="/Users/tomparee/Desktop/recombinationFraction.png", plot=p, height=1.6, width=6, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/recombinationFraction.pdf", plot=p,height=1.6, width=6, dpi=600)








load(file=paste0( "./Data/rec1EE/mixing_efficiency/predictedInterference_Ne1000_s0.01_h0.5_nloci", 1000, ".Rdata" ))
meanHR=aggregate(Er2_HillRobertson ~ nboot+rec+nsampledSNP,interference, mean)
meanHR$rec=factor(meanHR$rec, levels = c('wt', 'mut'), labels = c("Wild-type", "Mutant"))


nlowrecombing = aggregate(n~rec, subset(histrf, rf2<=0.01), sum)

plowrf = ggplot(nlowrecombing, aes(factor(rec, levels = c('wt', 'mut'), labels = c("Wild-type", "Mutant")), n))+
  theme_Publication3()+
  geom_bar(stat="identity", color="black", aes(fill=rec),width = 0.65/.pt,size=1/.pt)+
  scale_fill_manual(values=c(blueWT, orangeMut))+
  theme(legend.position = "none")+
  scale_y_continuous(expand=c(0,0), label=scientific_10)+
  xlab("")+ylab("Number of low \n recombining SNV")


pHR=ggplot(meanHR, aes(rec, Er2_HillRobertson, fill=rec), color='black')+
  stat_summary(fun.y = mean, 
               geom = "bar", size=0.7/.pt, color='black',width = 0.65/.pt)+ 
  stat_summary(fun.y = mean, 
               geom = "point", size=2/.pt)+ 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)*1.96, 
               fun.ymax = function(x) mean(x) + sd(x)*1.96, 
               geom = "errorbar",
               width=0.3/.pt,
               size=1/.pt,
               position = position_dodge(width=0.17))+
  scale_fill_manual(values=c(blueWT, orangeMut),name="")+
  theme_Publication3()+
  labs(y=expression(E(r^{2})))+
  xlab("")+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = "none")#+
#scale_y_continuous(expand = c(0.3,0))


p=grid.arrange(plowrf, pHR, widths = c(1.1,1))


ggsave(filename="/Users/tomparee/Desktop/Figures/experctedr2_2.pdf", plot=p, height=1.22, width=2.4, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/Figures/experctedr2_2.png", plot=p, height=1.22, width=2.4, dpi=600)



