#0.001560353
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(ggplot2)
library(gridExtra)

source("./basics_TP.R")

path = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/Data/rec1EE/Rpoly/"

load(file=paste0(path, 'Rpoly_all_rec1freq.Rdata'))
Rpoly = subset(Rpoly, !is.na(replicate))


scoef = do.call(rbind, lapply(split(Rpoly, Rpoly$replicate), function(x){
  s = try(summary( glm(pred_ratio~generation, data=x, family="binomial") )$coef[2,1])
  afc = try(summary(lm(pred_ratio~generation, data=x))$coef[2,1])
  if(class(s)[1]=='try-error'){s = NA}
  if(class(afc)[1]=='try-error'){afc = NA}
  x
  x = x[1, -which(colnames(x) %in% c("generation", "pred_ratio", "sample"))]
  x$scoef = as.numeric(s)
  x$afc =  as.numeric(afc)
  x
  
}))



#scoef$N2 = scoef$N
#scoef$N2[scoef$repro=='selfing'] = scoef$N2[scoef$repro=='selfing'] - 1000
#scoef$N2[scoef$repro=='outcrossing' & scoef$env=='salt' & scoef$N==10000] = scoef$N2[scoef$repro=='outcrossing' & scoef$env=='salt' & scoef$N==10000] + 1000



##################################################
### STATS #########################
library(emmeans)
library("lme4")




model1 = lm(scoef~as.factor(N)+env,  subset(scoef, repro=='outcrossing'))
model2 = lmer(scoef~as.factor(N)+env+(1|block),  subset(scoef, repro=='outcrossing'))
emm1 <- emmeans(model1, ~ N)
summary(emm,infer = c(TRUE, TRUE))

emm2 <- emmeans(model2, ~ N)
summary(emm2,infer = c(TRUE, TRUE))

# Model on s coef
model = lm(scoef~env + as.factor(N),  subset(scoef, repro=='outcrossing'))
summary(model,infer = c(TRUE, TRUE))

round(mean(subset(scoef, repro=='outcrossing' & N==1000)$scoef, na.rm=T), digits = 3)

model = lm(scoef~as.factor(N),  subset(scoef, repro=='outcrossing'))


mglm = glm(pred_ratio~ generation*as.factor(N),  subset(Rpoly, repro=='outcrossing'), family = "quasibinomial")
summary(mglm)


modelall = lmer(scoef~as.factor(N)+env+(1|block),  subset(scoef, repro=='outcrossing'))
modelnoenv = lmer(scoef~as.factor(N)+(1|block),  subset(scoef, repro=='outcrossing'))
modelnoN = lmer(scoef~env+(1|block),  subset(scoef, repro=='outcrossing'))
modelnointercept =  lmer(scoef~as.factor(N)+env+(1|block)-0,  subset(scoef, repro=='outcrossing'))

anova( lm(scoef~1,  subset(scoef, repro=='outcrossing' & N==1000)),
            lm(scoef~0,  subset(scoef, repro=='outcrossing' & N==1000)))



summary(modelall)
anova(modelall,modelnoenv)
anova(modelall,modelnoN)
anova(modelall,modelnointercept)

summary(model)
anova(model)

model = lm(scoef~ repro,  subset(scoef, env=='salt' & block == "Rpoly1"))
summary(model)


var.test(scoef~ repro,  subset(scoef, env=='salt' & block == "Rpoly1"),
         alternative = "two.sided")$p.value



var(subset(scoef, env=='salt' & block == "Rpoly1" & repro=='outcrossing')$scoef)
var(subset(scoef, env=='salt' & block == "Rpoly1" & repro=='selfing')$scoef)


# Model on generation
mglm = glm(pred_ratio~ generation*repro,  subset(Rpoly, env=='salt' & block == "Rpoly1"), family = "quasibinomial")
summary(mglm)


summary(glmer(pred_ratio~ generation*repro + (1|replicate),  subset(Rpoly, env=='salt' & block == "Rpoly1")))
mlmer = lmer(pred_ratio~ generation*repro + (1|replicate),  subset(Rpoly, env=='salt' & block == "Rpoly1"))
summary(mlmer)

int.outcross = summary(mlmer)$coef[1,1]
slp.outcross = summary(mlmer)$coef["generation",1]

slp.self = summary(mlmer)$coef["generation",1] + summary(mlmer)$coef["generation:reproselfing",1]
int.self = summary(mlmer)$coef[1,1] +  summary(mlmer)$coef["reproselfing",1]

Rpolyself.slopes = data.frame(intercept = c(int.outcross, int.self),
           slope = c(slp.outcross, slp.self),
           repro = c("outcross", "selfing"))

Rpolyself.slopes$x1 = 0
Rpolyself.slopes$x2 = 16
Rpolyself.slopes$y1 = Rpolyself.slopes$intercept + Rpolyself.slopes$slope * Rpolyself.slopes$x1
Rpolyself.slopes$y2 = Rpolyself.slopes$intercept + Rpolyself.slopes$slope * Rpolyself.slopes$x2



model0=glmer(pred_ratio~ 1 + (1|block/replicate),  subset(Rpoly, repro=="outcrossing"), family='binomial')
model1=glmer(pred_ratio~ generation + (1|block/replicate),  subset(Rpoly, repro=="outcrossing"),family='binomial')
anova(model0,model1,test="LRT")
summary(model1)

model0=lmer(pred_ratio~ 1 + (1|block/replicate),  subset(Rpoly, repro=="outcrossing"))
model1=lmer(pred_ratio~ generation + (1|block/replicate),  subset(Rpoly, repro=="outcrossing"))
anova(model0,model1,test="LRT")
summary(model1)

mlmer = lmer(pred_ratio~ generation*repro + (1|replicate),  subset(Rpoly, env=='salt' & block == "Rpoly1"))
summary(mlmer)

#"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"

ps = ggplot()+theme_Publication3()+
  geom_point(data = subset(scoef, repro=='outcrossing'), 
             aes(as.factor(N), scoef, color = as.factor(N)), position=position_jitter(width=0.1),size=1.75/.pt)+
  stat_summary(data = subset(scoef, repro=='outcrossing'), aes(as.factor(N), scoef),
               fun.y = mean,
               fun.ymin = function(x) mean(x) - 1.386*se(x), 
               fun.ymax = function(x) mean(x) + 1.386*se(x), 
               geom = "errorbar",linewidth=1.3/.pt, width = 0.1)+
  stat_summary(data = subset(scoef, repro=='outcrossing'), aes(as.factor(N), scoef),
               fun = mean,size=1/.pt)+
  facet_wrap(~factor(env, levels= c("standard", 'salt'), labels = c("Domestication", "Novel")), scales="free")+
  geom_hline(yintercept = 0, linetype='dashed', size=1/.pt)+
  xlab("Population size")+ylab(expression(bold("Selection of the "*bolditalic("rec-1")*" mutant")))+
  theme(axis.title = element_text(size=7,colour = "black", face = "bold"),    
        axis.text = element_text(size=7, colour = "black"),
        axis.text.x = element_text(size=7, colour = "black", face = "bold"),
        strip.text = element_text(size = 7, color = "black", face = "bold"),
        plot.title = element_text(size = 7, hjust = -0.05, face = "bold"),
        legend.position = "none")+
  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
  theme(legend.position = "none")+
  scale_x_discrete(breaks = c("1000","10000"), labels = c(bquote(bold(10^3)),bquote(bold(10^4))) )+
  ylim(-0.02,0.07)


mean(subset(scoef, repro=='outcrossing' & env=='High salt' & N==1000)$scoef,na.rm=T)-(se(subset(scoef, repro=='outcrossing' & env=='High salt' & N==1000)$scoef,na.rm=T)*1.386)

scoef$env = factor(scoef$env, levels= c("standard", 'salt'), labels = c("Domestication", "High salt"))
ps = ggplot()+theme_Publication3()+
  geom_point(data = subset(scoef, repro=='outcrossing'), 
             aes(env, scoef, color = as.factor(N)), position=position_jitter(width=0.1),size=1.2/.pt)+
  stat_summary(data = subset(scoef, repro=='outcrossing'), aes(env, scoef),
               fun.y = mean,
               fun.ymin = function(x) mean(x) - 1.96*se(x), 
               fun.ymax = function(x) mean(x) + 1.96*se(x), 
               geom = "errorbar",linewidth=1/.pt, width = 0.1)+
  stat_summary(data = subset(scoef, repro=='outcrossing'), aes(env, scoef),
               fun = mean,size=0.25/.pt)+
  facet_wrap(~N, scales="free")+
  geom_hline(yintercept = 0, linetype='dashed', size=1/.pt)+
  xlab("")+ylab(bquote(Indirect ~ selection ~ (s[i])))+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = seq(-0.1,0.1,0.01), lim=c(-0.02,0.07))




ggsave(filename="/Users/tomparee/Desktop/Figures/indirectselection.png", plot=ps, height=1.2, width=2.4, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/Figures/indirectselection.pdf", plot=ps, height=1.2, width=2.4, dpi=600)




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

# load( "/Users/tomparee/Documents/Documents - MacBook Pro de tom/Data/rec1EE/Rpoly/expectedLinkedSelection_estimatedFromCeMEE.Rdata")
# 
# meanafc = setNames(cbind(aggregate(afc~N, scoef, mean), aggregate(afc~N, scoef, se)[,"afc"]), c("N","afc","se"))
# meanafc$y = c(4,9)
# expectedLinkedSelection = subset(expectedLinkedSelection, boot != "obs" & winsize==1.5e6)
# pafc = ggplot()+theme_Publication3()+theme(legend.position = "none")+
#   geom_density(data=expectedLinkedSelection, aes(log(afc, base=10)), size = 1/.pt, fill="grey")+
#   scale_x_continuous(breaks = seq(-5,0,1), labels = expression(10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1), expand = c(0,0), lim=c(-5.1,-1.9))+
#   scale_y_continuous(expand = c(0,0), breaks = c(0,5,10))+
#   geom_point(data=meanafc, aes(log(afc, base = 10), y, color=as.factor(N)), size=3/.pt)+
#   geom_errorbarh(data=meanafc, aes(xmin = log(afc-(1.96*se), base=10), xmax=log(afc+(1.96*se), base=10), y=y, color=as.factor(N)), height = 1.5/.pt, size=1/.pt)+
#   scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
#   #xlab(expression(paste(Delta, bolditalic("p"))))+
#   xlab(expression(paste(bolditalic("rec-1"),bold(" mutant frequency change"))))+
#   ylab("Density")+
#   theme(axis.title = element_text(size=7,colour = "black", face = "bold"),    
#         axis.text = element_text(size=7, colour = "black"),
#         axis.text.x = element_text(size=7, colour = "black", face = "bold"),
#         strip.text = element_text(size = 7, color = "black", face = "bold"))
# 
# 
# p=grid.arrange(ps+theme(plot.margin = unit(c(5,5,35,5), "pt")), pafc, ncol=1, heights = c(2.9,1.1))
# 
# ggsave(filename="/Users/tomparee/Desktop/indirectselection.png", plot=p, height=2.9, width=3.2, dpi=600)
# ggsave(filename="/Users/tomparee/Desktop/indirectselection.pdf", plot=p, height=2.9, width=3.2, dpi=600)








# load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/Data/rec1EE/expectedLinkedSelection_estimatedFromRhom.Rdata")
# 
# meanscoef = setNames(cbind(aggregate(scoef~N+env, scoef, mean), aggregate(scoef~N+env, scoef, se)[,"scoef"]), c("N","env","scoef","se"))
# 
# meanscoef$y = c(0.5,1.5,1.5, 4.5)
# 
# 
# 
# expectedLinkedSelection = subset(expectedLinkedSelection, boot != "obs")
# #meanscoef=do.call(rbind, lapply(unique(expectedLinkedSelection$winsize), function(w){
# #  x=meanscoef
# #  x$winsize=w
# #  x
# #}))
# 
# 
# 
# 
# 
# psalt = ggplot(subset(expectedLinkedSelection, slink_salt>0))+theme_Publication3()+
#   geom_density(aes(log(slink_salt, base=10), group=winsize,fill=winsize), size = 1/.pt, alpha=0.85)+
#   geom_point(data=subset(meanscoef,env=='salt'), aes(log(scoef, base = 10), y, color=as.factor(N)), size=2/.pt)+
#   geom_errorbarh(data=subset(meanscoef,env=='salt'), 
#                  aes(xmin = log(scoef-(1.96*se), base=10), xmax=log(scoef+(1.96*se), base=10), y=y, color=as.factor(N)),
#                  height = 0.4/.pt, size=1/.pt)+
#   scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
#   scale_fill_gradient(low="grey85", high = "grey45")+
#   scale_x_continuous(breaks = seq(-7,0,1), labels = expression(10^-7,10^-6,10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1), expand = c(0,0), lim=c(-7.1,-1))+
#   scale_y_continuous(expand = c(0,0), breaks = seq(0,10,2))+
#   xlab(expression(paste(italic("rec-1")," mutant selection coeficient")))+
#   ylab("Density")+
#   ggtitle("High salt")
# 
# pdom = ggplot(subset(expectedLinkedSelection, slink_domestication>0))+theme_Publication3()+
#   geom_density(aes(log(slink_domestication, base=10), group=winsize,fill=winsize), size = 1/.pt, alpha=0.85)+
#   geom_point(data=subset(meanscoef, env=="standard"), aes(log(scoef, base = 10), y, color=as.factor(N)), size=2/.pt)+
#   geom_errorbarh(data=subset(meanscoef, env=="standard"), 
#                  aes(xmin = log(scoef-(1.96*se), base=10), xmax=log(scoef+(1.96*se), base=10), y=y, color=as.factor(N)),
#                  height = 1/.pt, size=1/.pt)+
#   scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
#   scale_fill_gradient(low="grey85", high = "grey45")+
#   scale_x_continuous(breaks = seq(-7,0,1), labels = expression(10^-7,10^-6,10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1), expand = c(0,0), lim=c(-7.1,-1))+
#   scale_y_continuous(expand = c(0,0), breaks = seq(0,10,2))+
#   xlab(expression(paste(italic("rec-1")," mutant selection coeficient")))+
#   ylab("Density")+
#   ggtitle("Domestication")
# 
# 
# p=gridExtra::grid.arrange(pdom+theme(legend.position="none"),psalt+theme(legend.position="none"),ncol=1)
# 
# ggsave(filename="/Users/tomparee/Desktop/Figures/expected_slink.png", plot=p, height=1.76, width=1.6, dpi=600)
# ggsave(filename="/Users/tomparee/Desktop/Figures/expected_slink.png.pdf", plot=p, height=1.76, width=1.6, dpi=600)
# 
# #ggplot(subset(expectedLinkedSelection, slink_salt>0))+theme_Publication3()+theme(legend.position = "none")+
# #  facet_grid(winsize~.)+
# #  geom_density(aes(log(slink_salt, base=10), group=winsize), size = 1/.pt, fill="grey")+
# #  geom_rect(data=meanscoef, aes(xmin = log(scoef-(1.96*se), base=10),
# #                                xmax= log(scoef+(1.96*se), base=10),
# #                                ymin=-Inf, ymax = Inf, fill=as.factor(N)), color=NA, alpha =0.1)+
# #  scale_fill_manual(values = c("#91D1C2FF","#00A087FF"))+
# #  scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
# #  geom_vline(data=meanscoef, aes(xintercept = log(scoef,base=10), color = as.factor(N)))+
# #  coord_cartesian(clip="off")
# 
#   
#   
#   
#   scale_x_continuous(breaks = seq(-5,0,1), labels = expression(10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1), expand = c(0,0), lim=c(-5.1,-1.9))+
#   scale_y_continuous(expand = c(0,0), breaks = c(0,5,10))+
#   geom_point(data=meanafc, aes(log(afc, base = 10), y, color=as.factor(N)), size=3/.pt)+
#   geom_errorbarh(data=meanafc, aes(xmin = log(afc-(1.96*se), base=10), xmax=log(afc+(1.96*se), base=10), y=y, color=as.factor(N)), height = 1.5/.pt, size=1/.pt)+
#   scale_color_manual(values = c("#91D1C2FF","#00A087FF"))+
#   #xlab(expression(paste(Delta, bolditalic("p"))))+
#   xlab(expression(paste(bolditalic("rec-1"),bold(" mutant frequency change"))))+
#   ylab("Density")+
#   theme(axis.title = element_text(size=7,colour = "black", face = "bold"),    
#         axis.text = element_text(size=7, colour = "black"),
#         axis.text.x = element_text(size=7, colour = "black", face = "bold"),
#         strip.text = element_text(size = 7, color = "black", face = "bold"))
# 
# 
# 
# 
# 
