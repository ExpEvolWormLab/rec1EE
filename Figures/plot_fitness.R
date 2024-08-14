PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")
library(emmeans)
library(data.table)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)

load(file = "adaptation_fitness/fitnessGain.Rdata") # load FITNESS (a list of data.frame)


##########################################################
############ Adaptation (fitness gain) ##############################

padaptsalt = ggplot()+
  theme_Publication3()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, color = "grey30", size=1/.pt)+
  facet_grid(.~factor(env,levels = c("domestication", "salt"), labels = c("Domestication", "High salt")),scales="free", space = "free")+
  theme(panel.margin.y=unit(20, "pt"))+
  geom_point(data=subset(FITNESS$fitness.gain, env=="salt"), aes(pop_random, fitness.gain.rel,color=rec), position = position_jitter(width = 0.25),  alpha = 0.2, size = 1.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.fitness.gain.pop,env=="salt"), aes(x=pop_random, ymin=lower.CL, ymax=upper.CL, color=rec),size = 1.3/.pt, width = 0.6/.pt)+
  geom_point(data=subset(FITNESS$mean.fitness.gain.pop,env=="salt"), aes(pop_random, fitness.gain.rel, color=rec), fill = "white",shape=21, size = 4/.pt,stroke = 2/.pt)+
  geom_point(data=subset(FITNESS$mean.fitness.gain.rec, env=='salt'), aes(xticks+0.25, fitness.gain.rel), size=5.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.fitness.gain.rec,env=="salt"), aes(x=xticks+0.25, ymin=lower.CL, ymax=upper.CL),size = 1.7/.pt, width = 2/.pt)+
  scale_color_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut))+
  scale_fill_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut))+
  scale_x_continuous(breaks = subset(FITNESS$mean.fitness.gain.rec, env=='salt')$xticks, labels = c("Wild-type", "Mutant"), 
                     limits = range(subset(FITNESS$mean.fitness.gain.pop,env=="salt")$pop_random)+c(-3,3))+
  ylim(range(FITNESS$fitness.gain$fitness.gain.rel,na.rm=T))+
  ylab("\u0394w (%)")



padaptngm = ggplot()+
  theme_Publication3()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, color = "grey30", size=1/.pt)+
  facet_grid(.~factor(env,levels = c("domestication", "salt"), labels = c("Domestication", "High salt")),scales="free", space = "free")+
  theme(panel.margin.y=unit(20, "pt"))+
  geom_point(data=subset(FITNESS$fitness.gain, env=="domestication"), aes(pop_random, fitness.gain.rel,color=rec), position = position_jitter(width = 0.25),  alpha = 0.2, size = 1.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.fitness.gain.pop,env=="domestication"), aes(x=pop_random, ymin=lower.CL, ymax=upper.CL, color=rec),size = 1.3/.pt, width = 0.6/.pt)+
  geom_point(data=subset(FITNESS$mean.fitness.gain.pop,env=="domestication"), aes(pop_random, fitness.gain.rel, color=rec), fill = "white",shape=21, size = 4/.pt,stroke = 2/.pt)+
  geom_point(data=subset(FITNESS$mean.fitness.gain.rec, env=="domestication"), aes(xticks+0.25, fitness.gain.rel), size=5.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.fitness.gain.rec,env=="domestication"), aes(x=xticks+0.25, ymin=lower.CL, ymax=upper.CL),size = 1.7/.pt, width = 2/.pt)+
  scale_color_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut))+
  scale_fill_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut))+
  scale_x_continuous(breaks = subset(FITNESS$mean.fitness.gain.rec, env=='domestication')$xticks, labels = c("Wild-type", "Mutant"), 
                     limits = range(subset(FITNESS$mean.fitness.gain.pop,env=="domestication")$pop_random)+c(-3,3))+
  ylim(range(FITNESS$fitness.gain$fitness.gain.rel,na.rm=T))+
  ylab("\u0394w (%)")




padaptV = grid.arrange(padaptngm,
                      padaptsalt,
                      nrow=2)



##############################################################
# Adaptation in salt without replicates dots (for HT) ###########

padaptsalt2 = ggplot()+
  theme_Publication3()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept = 0, color = "grey30", size=1/.pt)+
  theme(panel.margin.y=unit(20, "pt"))+
  geom_errorbar(data=subset(FITNESS$mean.fitness.gain.pop,env=="salt"), aes(x=pop_random, ymin=lower.CL, ymax=upper.CL, color=rec),size = 1.3/.pt, width = 1.2/.pt)+
  geom_point(data=subset(FITNESS$mean.fitness.gain.pop,env=="salt"), aes(pop_random, fitness.gain.rel, color=rec), fill = "white",shape=21, size = 4/.pt,stroke = 2/.pt)+
  geom_point(data=subset(FITNESS$mean.fitness.gain.rec, env=='salt'), aes(xticks+0.25, fitness.gain.rel), size=5.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.fitness.gain.rec,env=="salt"), aes(x=xticks+0.25, ymin=lower.CL, ymax=upper.CL),size = 1.7/.pt, width = 2/.pt)+
  scale_color_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut))+
  scale_fill_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut))+
  scale_x_continuous(breaks = subset(FITNESS$mean.fitness.gain.rec, env=='salt')$xticks, labels = c("Wild-type", "Mutant"), 
                     limits = range(subset(FITNESS$mean.fitness.gain.pop,env=="salt")$pop_random)+c(-3,3))+
  ylab("\u0394w (%)")


ggsave(filename="Figures/fitness_HT.png", plot=padaptsalt2, height=1.1, width=1.7, dpi=600)
ggsave(filename="Figures/fitness_HT.pdf", plot=padaptsalt2, height=1.1, width=1.7, dpi=600)



#pancestral = ggplot(data=FITNESS$mean.ancestral.fitness)+
#  theme_Publication3()+theme(legend.position = "none")+
#  theme(axis.title.x = element_blank())+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#  facet_grid(.~factor(env,levels = c("domestication", "salt"), labels = c("Domestication", "High salt")))+
#  theme(panel.margin.y=unit(20, "pt"))+
#  geom_point(data=FITNESS$ancestal.fitness, aes(rec,rel.fitness, color=rec), position = position_jitter(width = 0.1), alpha = 0.25, size = 1.5/.pt)+
#  geom_errorbar(data=FITNESS$mean.ancestral.fitness, aes(x=rec, ymin=lower.CL, ymax=upper.CL, color=rec),size = 1.3/.pt, width = 0.6/.pt)+
#  geom_point(data=FITNESS$mean.ancestral.fitness, aes(rec, rel.fitness,color=rec), fill = "white",shape=21, size = 4/.pt,stroke = 3/.pt)+
#  geom_segment(data=subset(FITNESS$mean.fitness.gain.pop,env=="salt"), x=-Inf, xend = -Inf, y=-Inf, yend= Inf, size = 1/.pt)+
#  scale_color_manual(values = c(blueWT, orangeMut))+
#  scale_x_discrete(breaks = c("wt", "mut"), labels = c("Wild-type", "Mutant"))+
#  ylab("Ancestral fitness")


##########################################################
############ Ancestral gain ##############################

pancestralngm = ggplot(data=subset(FITNESS$mean.ancestral.fitness, env=="domestication"))+
  theme_Publication3()+theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~factor(env,levels = c("domestication", "salt"), labels = c("Domestication", "High salt")), nrow=2)+
  theme(panel.margin.y=unit(20, "pt"))+
  geom_point(data=subset(FITNESS$ancetral.fitness, env=="domestication"), aes(rec,rel.fitness, color=rec), position = position_jitter(width = 0.1), alpha = 0.25, size = 1.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.ancestral.fitness, env=="domestication"), aes(x=rec, ymin=lower.CL, ymax=upper.CL, color=rec),size = 1.3/.pt, width = 0.6/.pt)+
  geom_point(data=subset(FITNESS$mean.ancestral.fitness, env=="domestication"), aes(rec, rel.fitness,color=rec), fill = "white",shape=21, size = 4/.pt,stroke = 3/.pt)+
  scale_color_manual(values = c(blueWT, orangeMut))+
  scale_x_discrete(breaks = c("wt", "mut"), labels = c("Wild-type", "Mutant"))+
  ylab("Ancestral fitness (w)")+
  ylim(range(FITNESS$ancetral.fitness$rel.fitness,na.rm = T))

pancestralsalt = ggplot(data=subset(FITNESS$mean.ancestral.fitness, env=='salt'))+
  theme_Publication3()+theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~factor(env,levels = c("domestication", "salt"), labels = c("Domestication", "High salt")), nrow=2)+
  theme(panel.margin.y=unit(20, "pt"))+
  geom_point(data=subset(FITNESS$ancetral.fitness, env=='salt'), aes(rec,rel.fitness, color=rec), position = position_jitter(width = 0.1), alpha = 0.25, size = 1.5/.pt)+
  geom_errorbar(data=subset(FITNESS$mean.ancestral.fitness, env=='salt'), aes(x=rec, ymin=lower.CL, ymax=upper.CL, color=rec),size = 1.3/.pt, width = 0.6/.pt)+
  geom_point(data=subset(FITNESS$mean.ancestral.fitness, env=='salt'), aes(rec, rel.fitness,color=rec), fill = "white",shape=21, size = 4/.pt,stroke = 3/.pt)+
  scale_color_manual(values = c(blueWT, orangeMut))+
  scale_x_discrete(breaks = c("wt", "mut"), labels = c("Wild-type", "Mutant"))+
  ylab("Ancestral fitness (w)")+
  ylim(range(FITNESS$ancetral.fitness$rel.fitness,na.rm = T))

pancestralV = grid.arrange(pancestralngm, pancestralsalt,nrow=2)

#ggsave(filename="/Users/tomparee/Desktop/manhattan.png", plot=p, height=2.86, width=3.5, dpi=600)



p = grid.arrange(pancestral, padapt, widths = c(1,1.5))
ggsave(filename="/Users/tomparee/Desktop/fitness2.png", plot=p, height=1.15, width=5.5, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/fitness2.pdf", plot=p, height=1.15, width=5.5, dpi=600)

pV = grid.arrange(pancestralV, padaptV, widths = c(1,1.4))
ggsave(filename="/Users/tomparee/Desktop/Figures/fitnessV.png", plot=pV, height=2.45, width=3.12, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/Figures/fitnessV.pdf", plot=pV, height=2.45, width=3.12, dpi=600)



grid.arrange(pancestral, padapted, nrow = 1, widths=c(1,1.5))





