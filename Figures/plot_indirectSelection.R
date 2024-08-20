PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

library(ggplot2)
library(ggplot2)
library(gridExtra)

load("indirect_selection/MR_selectionCoeficients.Rdata")

#model <- lm(scoef~as.factor(N)+env, data = scoef)
#value = as.data.frame(emmeans(model, specs = ~ N+env, level= 0.95))



ps=ggplot()+theme_Publication3()+
  geom_point(data = scoef, 
             aes(env, scoef, color = as.factor(N)), position=position_jitter(width=0.1),size=1.2/.pt)+
  stat_summary(data = scoef, aes(env, scoef),
               fun.y = mean,
               fun.ymin = function(x) mean(x) - 1.96*se(x), 
               fun.ymax = function(x) mean(x) + 1.96*se(x), 
               geom = "errorbar",linewidth=1/.pt, width = 0.1)+
  stat_summary(data = scoef, aes(env, scoef),
               fun = mean,size=0.25/.pt)+
  facet_wrap(~factor(N, levels = c(10000,1000)), scales="free")+
  geom_hline(yintercept = 0, linetype='dashed', size=1/.pt)+
  xlab("")+ylab(bquote(Indirect ~ selection ~ (s[i])))+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = seq(-0.1,0.1,0.02), lim=c(-0.02,0.07))



ggsave(filename="Figures/indirectSelection.png", plot=ps, height=1.4, width=2.7, dpi=600)
ggsave(filename="Figures/indirectSelection.pdf", plot=ps,  height=1.4, width=2.7, dpi=600)


summary(lm(scoef~1, data=subset(scoef, N==10000))) # 0.00041 ***
summary(lm(scoef~1, data=subset(scoef, N==10000 & env=="domestication"))) #0.00641 **
summary(lm(scoef~1, data=subset(scoef, N==10000 & env=="salt"))) # 0.031 *

summary(lm(scoef~1, data=subset(scoef, N==1000))) # 2.78e-07 ***
summary(lm(scoef~1, data=subset(scoef, N==1000 & env=="domestication"))) # 6.44e-05 ***
summary(lm(scoef~1, data=subset(scoef, N==1000 & env=="salt"))) #  0.00095 ***



# 
#ps2=ggplot(data = subset(scoef, env=="salt"), 
#       aes(as.factor(N), scoef, color = as.factor(N)))+theme_Publication3()+
#  geom_point(position=position_jitter(width=0.1),size=1.2/.pt)+
#  stat_summary(fun.y = mean,
#               fun.ymin = function(x) mean(x) - 1.96*se(x), 
#               fun.ymax = function(x) mean(x) + 1.96*se(x), 
#               geom = "errorbar",linewidth=1/.pt, width = 0.1, color='black')+
# stat_summary(fun = mean,size=0.25/.pt, color='black')+
#  geom_hline(yintercept = 0, linetype='dashed', size=1/.pt)+
#  xlab("Population size (N)")+ylab(bquote(Indirect ~ selection ~ (s[i])))+
#  theme(legend.position = "none")+
#  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
#  theme(legend.position = "none")+
#  scale_y_continuous(breaks = seq(-0.1,0.1,0.02), lim=c(-0.02,0.07))




load(file = "./Data/rec1EE/Rpoly/selectionEEV1403.Rdata")

pf = ggplot(freqEEV1403, aes(generation, freq))+
  geom_line(aes(group=marker), color='orange')+stat_summary(geom="line")+
  ylim(0,0.33)+facet_wrap(~env ,nrow=2, scales = 'free')+
  coord_cartesian(expand=0)+
  theme_Publication3()+
  ylab("EEV1403 alleles frequency")+xlab("Generation")



selectionEEV1403$env[selectionEEV1403$env=="dom"]='standard'
pss = ggplot()+theme_Publication3()+
  geom_point(data = selectionEEV1403, aes(x=1, y=s), color='orange', size=2.2/.pt)+
  geom_errorbar(data=selectionEEV1403, aes(x=1, ymin = s-(1.386*se), ymax = s+(1.386*se)), color='orange', width = 0.1)+
  geom_rect(data=selectionEEV1403, aes(xmin=-Inf, xmax=Inf,, ymin = s-(1.386*se), ymax = s+(1.386*se)), fill='orange', alpha=0.1, color=NA)+
  stat_summary(data = subset(scoef, repro=='outcrossing' & N == 1000), aes(1.5, scoef),
               fun.y = mean,
               fun.ymin = function(x) mean(x) - 1.386*se(x), 
               fun.ymax = function(x) mean(x) + 1.386*se(x), 
               geom = "errorbar", width = 0.1, color="#91D1C2FF")+
  stat_summary(data = subset(scoef, repro=='outcrossing' &  N == 1000), aes(1.5, scoef),
               fun = mean,size=0.25/.pt, color="#91D1C2FF")+
  stat_summary(data = subset(scoef, repro=='outcrossing' & N == 10000), aes(2, scoef),
               fun.y = mean,
               fun.ymin = function(x) mean(x) - 1.386*se(x), 
               fun.ymax = function(x) mean(x) + 1.386*se(x), 
               geom = "errorbar", width = 0.1,color = "#00A087FF")+
  stat_summary(data = subset(scoef, repro=='outcrossing' &  N == 10000), aes(2, scoef),
               fun = mean,size=0.25/.pt, color = "#00A087FF")+
  facet_wrap(~factor(env, levels = c('standard','salt')) ,nrow=2, scales = 'free')+
  coord_cartesian(ylim=c(0,0.03))+
  ylab("Selection coeficient")+
  xlab("")+
  theme(axis.ticks.x = element_line(color='white'), axis.text.x = element_text(color='white'))

p=grid.arrange(pf, pss, ncol=2)

ggsave(filename="/Users/tomparee/Desktop/Figures/selectionEEV1403.png", plot=p, height=2.57, width=2.68, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/Figures/selectionEEV1403.pdf", plot=p, height=2.57, width=2.68, dpi=600)



