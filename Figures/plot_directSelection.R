
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

library(gridExtra)
library(ggplot2)

load("direct_selection/direct_selection_bootstrap.Rdata")
load("direct_selection/direct_selection_rec1freq.Rdata")


## determine p-value (different from 0) given the overlap between 0 and the bootsrap distribution
#pvals = apply(CI[[2]], 2, function(x){ (2*min(sum(x<0),sum(x>0)))/sum(!is.na(x))})

#scoef = do.call(rbind, lapply(split(compet, compet$id), function(x){
#  cbind( unique(x[,c("id", "env", "genotype", "initialfreq_aimed", "malefreq", "block", "repro")]), 
#         data.frame(s=glm(rec1freq~generation, x, family="quasibinomial")$coef[2]))
#}))
#save(scoef, file="direct_selection/direct_selection_scoef.Rdata")

load("direct_selection/direct_selection_scoef.Rdata")

pfreq=ggplot(subset(compet, initialfreq_aimed==0.5), aes(generation, rec1freq, color = genotype))+
  theme_Publication3()+
  theme(legend.position = "none")+
  stat_summary(aes(group = sample),geom="line", alpha=0.7, size=1/.pt)+
  facet_wrap(~factor(env, levels = c("ngm", "nacl"), labels = c("Domestication", "High salt")),nrow=1, scales="free")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(lim=c(0.35,0.7))+
  scale_color_manual(values = c("#4DBBD5FF","#3C5488FF"))+
  ylab(expression(paste(italic("rec-1"), "mutant frequency")))+
  xlab("Generation")




pdis=ggplot()+theme_Publication3()+
    geom_density(data=scoef, aes(s), fill="grey90", size=1/.pt)+
    geom_vline(xintercept = 0, linetype='dashed', size=1/.pt)+
    geom_point(data = subset(CI[[1]], test == "gen"),
               aes(x=s, y=3), size=1.3/.pt)+
    geom_errorbarh(data = subset(CI[[1]], test == "gen"), aes(xmin=low.CI95, xmax=high.CI95, y=3),height=0.45, size=1.4/.pt)+
    coord_cartesian(expand = 0)+
    xlab(bquote(Direct ~ selection ~ (s[d]))) +
    ylab('Density')+theme(strip.text = element_text(size=7))

# just change some names
#CI[[1]]$test2= factor(CI[[1]]$test, levels = c("gen", "env", "genotype", "repro"), c("Global","Environment", "Genotype", "Sex ratio"))
#posy = c(4, 0.7,1.5,2.2)
#pdis=ggplot()+theme_Publication3()+
#  geom_density(data=scoef, aes(s), fill="grey90", size=1/.pt)+
#  geom_vline(xintercept = 0, linetype='dashed', size=1/.pt)+
#  geom_point(data = subset(CI[[1]], test %in% c("gen","repro", "genotype", "env")),
#             aes(x=s, y=posy, color=test2, size=test2))+
#  scale_colour_manual(breaks =  c("Global","Environment", "Genotype", "Sex ratio"), values = c("black", cbs[c(3,1,2)]), name="")+
#  scale_size_manual(breaks =  c("Global","Environment", "Genotype", "Sex ratio"), values = c(1.3/.pt, rep(0.6/.pt,3)), name="")+
#  geom_errorbarh(data = subset(CI[[1]], test %in% c("gen","repro", "genotype", "env")), aes(xmin=low.CI95, xmax=high.CI95, y=posy, color=test2, size=test2),height=0.3, size=0.8/.pt)+
#  geom_errorbarh(data = subset(CI[[1]], test == "gen"), aes(xmin=low.CI95, xmax=high.CI95, y=posy[1]),height=0.45, size=1.4/.pt)+
#  coord_cartesian(expand = 0)+
#  xlab(bquote(Direct ~ selection ~ (s[d]))) +
#  ylab('Density')+theme(strip.text = element_text(size=7))

p=grid.arrange(pfreq+theme(plot.margin=unit(c(5,10,5,20), "pt")),
               pdis+theme(plot.margin=unit(c(20,5,5,10), "pt"))+theme(legend.position = "none"), nrow=1, widths = c(1.5,1))

ggsave(filename="Figures/directselection.png", plot=p, height=1.2, width=3.8, dpi=600)
ggsave(filename="Figures/directselection.pdf", plot=p, height=1.2, width=3.8, dpi=600)


effects = CI[[1]]
effects = effects[match(c("EEV1401","EEV1402", "ngm","nacl", "selfing","outcrossing","gen"),effects$test), ]
effects$test2 = c("EEV1401","EEV1402", "Domestication","High salt", "Low","High","")
effects$group = c("Genetic background", "Genetic background", "Environment","Environment", "Male frequency", "Male frequency", "Average")
peff = ggplot(effects)+theme_Publication3()+theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype='dashed', linewidth = 1/.pt)+
  geom_point(aes(test2, s,color=test), size=1.6/.pt)+
  geom_errorbar(aes(x=test2,ymin=low.CI95,ymax=high.CI95, color=test), width=0.5/.pt, size=1/.pt)+
  facet_grid(~group, scale='free_x', space='free')+
  scale_color_manual(breaks = c("EEV1401","EEV1402", "ngm","nacl", "selfing","outcrossing","gen"),
                     values = c("#4DBBD5FF","#3C5488FF","#91D1C2FF", "#00A087FF", "#E69900","#CC6600", "black"))+
  ylab(bquote(Direct ~ selection ~ (s[d])))+xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggsave(filename="Figures/directselection_effects.png", plot=peff, height=1.5, width=3.1, dpi=600)
ggsave(filename="Figures/directselection_effects.pdf", plot=peff, height=1.5, width=3.1, dpi=600)


