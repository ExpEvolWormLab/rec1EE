library(readxl)
library(ggplot2)
library(data.table)
library(lme4)
library(lmerTest)
library(ggpubr)
library(vipor)
library(ggbeeswarm)
source("./basics_TP.R")
source('./Data/Experimental_Evolution_R/EE_Rhom/poolseq/utils_poolseq_analysis.R')

logit = qlogis
inv.logit=plogis


fitnessassay <- read_excel("./Data/rec1EE/fitness/fitnessassay_all.xlsx")





# 4 different experimental blocks:
# Different thaw and maintenance
# The GFPL1 tester is a differet thaw between blocks, but the same within blocks
# Blocks assayed on different days, except block C and D the same same day
fitnessassay$block[fitnessassay$block==1]="A"
fitnessassay$block[fitnessassay$block==2]="B"
fitnessassay$block[fitnessassay$block==3]="C"
fitnessassay$block[fitnessassay$block==4]="D"

fitnessassay$environment[fitnessassay$environment=="salt"]="nacl"

# Some identifier
fitnessassay$popexp=paste0(fitnessassay$populations,'_', fitnessassay$block,'_',fitnessassay$code)
fitnessassay$rec = factor(fitnessassay$rec, levels=c('wt','mut'))
fitnessassay$environment = factor(fitnessassay$environment, levels=c('ngm','nacl'))
fitnessassay$rec_evolution = paste0(fitnessassay$rec,fitnessassay$evolution)
fitnessassay$rec_evolution = factor(fitnessassay$rec_evolution, levels=c("wtancestor","wtevolved","mutancestor","mutevolved"))

# Initial mix GFP frequency 
fitnessassay$initial.f.gfp.l1 = (fitnessassay$initial.gfp.l1)/(fitnessassay$initial.all.l1)

# Selfing rate = 1 - 2*male_frequency
fitnessassay$self.rate = 1- 2*(fitnessassay$f.male)
fitnessassay$self.rate[fitnessassay$self.rate<0]=0 # minimum = 0

# GFP frequency after one generation competition
fitnessassay$f.gfp.l1 = (fitnessassay$gfp.l1)/(fitnessassay$all.l1)

# Initial male frequency
fitnessassay$initial.f.male=(fitnessassay$initial.male.tested)/(fitnessassay$initial.male.tested+fitnessassay$initial.herm.tested)
fitnessassay$initial.f.male[fitnessassay$initial.f.male>0.5]=0.5

# Take out populations SL3 and SHR4
# They have a significant drop in frequency during evolution
fitnessassay = subset(fitnessassay, !(populations %in% c("SLR3", "SHR4")))



# ggplot(data=fitnessassay)+stat_summary(aes(rec_evolution,logit(f.gfp.l1), group=popexp, color=rec), position = position_dodge(0.3))+stat_summary(aes(rec_evolution,logit(f.gfp.l1)))+facet_grid(block~environment,scales="free")
# 
# 
# # Look at distribution of GFP initial frequency
# ggplot()+geom_density(data=fitnessassay[fitnessassay$rec=='wt'& fitnessassay$tech.repl==1,], aes(logit(initial.f.gfp.l1)), color='black')+
#   geom_density(data=fitnessassay[fitnessassay$rec=='mut'& fitnessassay$tech.repl==1,], aes(logit(initial.f.gfp.l1)), color='red')+
#   facet_grid(environment~block,scales="free")
# 
# ggplot()+geom_density(data=fitnessassay[fitnessassay$rec=='wt'& fitnessassay$tech.repl==1,], aes(logit(initial.f.gfp.l1)), color='black')+
#   geom_density(data=fitnessassay[fitnessassay$rec=='mut'& fitnessassay$tech.repl==1,], aes(logit(initial.f.gfp.l1)), color='red')+
#   facet_grid(environment~.,scales="free")
# 
# aggregate(logit(initial.f.gfp.l1)~block, data=fitnessassay, mean)
# 
# # Look at the difference in initial GFP frequency on a logit scale, correcting for block effect
# # only to look, to used
# fitnessassay$initial.f.gfp.l1.norm = logit(fitnessassay$initial.f.gfp.l1)
# for(i in c("A","B","C","D")){
#   mix=mean(fitnessassay$initial.f.gfp.l1.norm[fitnessassay$block==i],na.rm=T)
#   fitnessassay$initial.f.gfp.l1.norm[fitnessassay$block==i]=fitnessassay$initial.f.gfp.l1.norm[fitnessassay$block==i]-mix
# }
# 
# ggplot()+geom_density(data=fitnessassay[fitnessassay$rec=='wt',], aes(initial.f.gfp.l1.norm), color='black')+
#   geom_density(data=fitnessassay[fitnessassay$rec=='mut',], aes(initial.f.gfp.l1.norm), color='red')+
#   facet_grid(environment~.)



##############################################################
#### Mating table to find the heterozigosity of GFP

fitnessassay$theo.het.gfp = NA
fitnessassay$theo.f.male.nonGFP = NA
fitnessassay$theo.f.male.GFP = NA
fitnessassay$theo.het.gfp.herm=NA
fitnessassay$rel.fitness=NA

for(i in 1:nrow(fitnessassay)){
  target = fitnessassay[i,]
  #F = S = selfing rate
  S = target$self.rate
  
  f0.nonGFP = 1 - fitnessassay$initial.f.gfp.l1[i]
  f0.GFP = fitnessassay$initial.f.gfp.l1[i]
  
  #p0.male.nonGFP = fitnessassay$initial.male.tested[i]
  #p0.herm.nonGFP = fitnessassay$initial.herm.tested[i]
  
  #if(p0.male.nonGFP>p0.herm.nonGFP) p0.male.nonGFP=p0.herm.nonGFP #max male freq = 0.5
  
  p0.male.GFP = fitnessassay$initial.male.gfp[i]
  p0.herm.GFP = fitnessassay$initial.herm.gfp[i]
  
  
  f0.male.nonGFP = fitnessassay$initial.f.male[i]
  f0.herm.nonGFP = 1-f0.male.nonGFP
  
  f0.male.GFP = p0.male.GFP/(p0.herm.GFP+p0.male.GFP)
  f0.herm.GFP = 1-f0.male.GFP
  
  
  
  #f2x = WT allele frequency in herm = f.WT * f.WT.h / (f.WT * f.WT.h + GFP.WT * GFP.WT.h)
  f2x = (f0.nonGFP*f0.herm.nonGFP)/((f0.nonGFP*f0.herm.nonGFP)+(f0.GFP*f0.herm.GFP))
  
  f2xx = (f0.nonGFP*f0.herm.nonGFP)/((f0.nonGFP*f0.herm.nonGFP)+(f0.GFP*f0.herm.GFP))
  
  f1x = 1 - f2x #herm
  f1xx = 1 - f2xx #male
  
  #Fitness in nonGFP (W2)
  
  a2 = f2x*f2xx*(1-S)
  b2 = f2x*S
  c2 = -1*(1 - target$f.gfp.l1)
  
  W2 = cbind((-b2 + sqrt(b2^2 - 4*a2*c2)) / (2*a2),
             (-b2 - sqrt(b2^2 - 4*a2*c2)) / (2*a2))[which.max(cbind((-b2 + sqrt(b2^2 - 4*a2*c2)) / (2*a2),
                                                                    (-b2 - sqrt(b2^2 - 4*a2*c2)) / (2*a2)))]
  
  p2x = f2x*W2
  p2xx = f2xx*W2
  
  #Fitness in GFP (W1)
  
  a1 = f1x*f1xx*(1-S)
  b1 = p2x*f1xx*(1-S) + p2xx*f1x*(1-S) + f1x*S
  c1 = p2x*p2xx*(1-S) + p2x*S - 1
  
  
  W1 = cbind((-b1 - sqrt(b1^2 - 4*a1*c1)) / (2*a1),
             (-b1 + sqrt(b1^2 - 4*a1*c1)) / (2*a1))[which.max(cbind((-b1 - sqrt(b1^2 - 4*a1*c1)) / (2*a1),(-b1 + sqrt(b1^2 - 4*a1*c1)) / (2*a1)))]
  
  #verify = f2xx*f2x*W2*W2*(1-S)+ W2*f2x*S + f1xx*f1x*W1*W1*(1-S)+ W1*f1x*S + f1x*f2xx*W1*W2*(1-S) + f1xx*f2x*W1*W2*(1-S)
  #should equal 1
  
  #rel.fitness = log((f0.nonGFP*W2)/(f0.GFP*W1)) - log(f0.nonGFP/f0.GFP)
  
  theo.het.gfp = (f1x*f2xx*W1*W2*(1-S) + f1xx*f2x*W1*W2*(1-S))/(f1xx*f1x*W1*W1*(1-S)+ W1*f1x*S + f1x*f2xx*W1*W2*(1-S) + f1xx*f2x*W1*W2*(1-S))
  theo.f.male.nonGFP = (0.5*(f2xx*f2x*W2*W2*(1-S)))/(f2xx*f2x*W2*W2*(1-S)+ W2*f2x*S)
  theo.f.male.GFP = (0.5*(f1xx*f1x*W1*W1*(1-S) + f1x*f2xx*W1*W2*(1-S) + f1xx*f2x*W1*W2*(1-S)))/(f1xx*f1x*W1*W1*(1-S)+ W1*f1x*S + f1x*f2xx*W1*W2*(1-S) + f1xx*f2x*W1*W2*(1-S))
  
  theo.het.gfp.herm= (0.5*0.5*(f1x*f2xx*W1*W2*(1-S)) + 0.5*0.5*(f1xx*f2x*W1*W2*(1-S)))/((W1*f1x*S)+(0.5*(f1x*f2xx*W1*W2*(1-S)) + 0.5*(f1xx*f2x*W1*W2*(1-S))))
  
  
  if(is.na(f0.nonGFP)|is.na(f0.male.nonGFP)|is.na(f0.herm.nonGFP)|is.na(f0.GFP)|is.na(f0.male.GFP)|is.na(f0.herm.GFP)|is.na(S)|is.na(target$f.gfp.l1)){
    theo.het.gfp = NA
    theo.f.male.nonGFP = NA
    theo.f.male.GFP = NA
    theo.het.gfp.herm=NA
    rel.fitness=NA
  }
  
  fitnessassay$theo.het.gfp[i] = theo.het.gfp
  fitnessassay$theo.het.gfp.herm[i]= theo.het.gfp.herm
  fitnessassay$theo.f.male.nonGFP[i] = theo.f.male.nonGFP
  fitnessassay$theo.f.male.GFP[i] = theo.f.male.GFP
  #fitnessassay$rel.fitness[i]=rel.fitness
  #theo = rbind(theo, cbind(W2,W1,theo.het.gfp,theo.f.male.nonGFP, theo.f.male.GFP, verify))
  
}

# Find the frequency of GFP haplotypes (f.chrom.gfp)
fitnessassay$factor.het.gfp = (((2*(1-fitnessassay$theo.het.gfp))+fitnessassay$theo.het.gfp)/2)
fitnessassay$f.chrom.gfp = fitnessassay$factor.het.gfp * (fitnessassay$gfp.l1/fitnessassay$all.l1)

# Rel fitness in ln(1-p1 / p1) - ln(1-p0 / p0), where p in the frequency of gfp haplotype
fitnessassay$rel.fitness = 1 + log((1-fitnessassay$f.chrom.gfp)/fitnessassay$f.chrom.gfp) - log((1-fitnessassay$initial.f.gfp.l1)/fitnessassay$initial.f.gfp.l1)



# Mean fitness in ancestors
ancestors.fitness = aggregate(rel.fitness~block+rec+environment, data=subset(fitnessassay, evolution == "ancestor"), mean)

# Calculate the fitness gain = rel fitness - respectibe ancestral rel fitness
fitnessassay$fitness.gain = NA
fitnessassay$fitness.gain.rel = NA
for(i in 1:nrow(ancestors.fitness)){
  env = ancestors.fitness$environment[i]
  block = ancestors.fitness$block[i]
  rec= ancestors.fitness$rec[i]
  fit = ancestors.fitness$rel.fitness[i]
  
  tar = which(fitnessassay$environment==env & fitnessassay$block==block & fitnessassay$rec==rec)
  
  fitnessassay$fitness.gain[tar]=fitnessassay$rel.fitness[tar] - fit
  fitnessassay$fitness.gain.rel[tar]=(fitnessassay$rel.fitness[tar]/fit)-1
}


fitnessassay = fitnessassay[,c("block","populations", "rec","environment", "evolution", "rel.fitness", "fitness.gain","fitness.gain.rel")]    
colnames(fitnessassay)= c("block","pop", "rec","env", "generation", "rel.fitness", "fitness.gain","fitness.gain.rel")
fitnessassay$generation[fitnessassay$generation=="evolved"] = "G40"
fitnessassay$generation[fitnessassay$generation=="ancestor"] = "G0"
fitnessassay$env=as.character(fitnessassay$env)
fitnessassay$env[which(fitnessassay$env=="nacl")] = "salt"

write.table(as.data.frame(fitnessassay), file = "./Data/rec1EE/Fitness/fitness.txt", row.names = F)


