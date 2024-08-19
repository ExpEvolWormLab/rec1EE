library(readr)
library(emmeans)
library(ggplot2)
library(gridExtra)
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")


############################
#### FERTILITY #############
saltfertility = read_csv("salt/Fertility_oocytes_salt.csv")


# Test fertility ~ nacl
modelfertility = glm(fertility~as.factor(nacl_concentration), data=saltfertility, family = "quasipoisson")
emmeans(modelfertility, pairwise ~ nacl_concentration)

#############################
#### MALE SUCCESS ###########

malefreq = read_csv(file="salt/f_male_salt.csv")


modelfmale = glm(cbind(nmale, nherm)~as.factor(nacl_concentration), data=malefreq, family = "binomial")
emmeans(modelfmale, pairwise ~ nacl_concentration)

mean_malefreq = as.data.frame(emmeans(modelfmale, ~ nacl_concentration, type="response"))


pm = ggplot(mean_malefreq)+
  geom_errorbar(size=1/.pt, width=0.3/.pt, aes(x=as.factor(nacl_concentration), ymin = asymp.LCL, ymax=asymp.UCL))+
  geom_point(size=1.75/.pt, aes(as.factor(nacl_concentration), prob))+theme_Publication3()+
  xlab('NaCl concentration (mM)')+ylab('Male frequency')+
  labs(color='Initial male\n frequency')+
  theme(legend.position="top")+
  theme(panel.border=element_rect(linetype=1, fill=NA))



pf = ggplot(saltfertility, aes(as.factor(nacl_concentration), fertility))+
  geom_boxplot(fill='grey90', outlier.size =1/.pt, size=1/.pt, width=1.5/.pt)+theme_minimal()+theme_Publication3()+
  xlab('NaCl concentration (mM)')+ylab('Hermaphrodite fertility')+
  theme(panel.border=element_rect(linetype=1, fill=NA))

p=grid.arrange(pm, pf, nrow=1, widths = c(1,1))

ggsave(filename="/Figures/salt_effect.png", plot=p, height=1.26, width=3, dpi=600)
ggsave(filename="/Figures/salt_effect.pdf", plot=p, height=1.26, width=3, dpi=600)

