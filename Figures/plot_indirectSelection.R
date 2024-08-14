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
  facet_wrap(~N, scales="free")+
  geom_hline(yintercept = 0, linetype='dashed', size=1/.pt)+
  xlab("")+ylab(bquote(Indirect ~ selection ~ (s[i])))+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = seq(-0.1,0.1,0.01), lim=c(-0.02,0.07))


# 
ps2=ggplot(data = subset(scoef, env=="salt"), 
       aes(as.factor(N), scoef, color = as.factor(N)))+theme_Publication3()+
  geom_point(position=position_jitter(width=0.1),size=1.2/.pt)+
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - 1.96*se(x), 
               fun.ymax = function(x) mean(x) + 1.96*se(x), 
               geom = "errorbar",linewidth=1/.pt, width = 0.1, color='black')+
  stat_summary(fun = mean,size=0.25/.pt, color='black')+
  geom_hline(yintercept = 0, linetype='dashed', size=1/.pt)+
  xlab("Census polunation size (N)")+ylab(bquote(Indirect ~ selection ~ (s[i])))+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = seq(-0.1,0.1,0.02), lim=c(-0.02,0.07))


ggsave(filename="Figures/indirectSelection_HT.png", plot=ps2, height=1.2, width=1.3, dpi=600)
ggsave(filename="Figures/indirectSelection_HT.pdf", plot=ps2, height=1.2, width=1.3, dpi=600)

