
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

load("heritability/h2~deltaR_0.9.Rdata")

h2 = subset(h2, region != "all")
p=ggplot(data=h2)+
  geom_bar(aes(region, h2, group = region), stat="identity", fill="lightgrey",color="black",width=0.65/.pt, size=1/.pt)+
  geom_point(aes(region, h2, group = region), size=1.5/.pt)+
  geom_errorbar(aes(x=region, ymin=h2-SE, ymax=h2+SE, group = region), width = 0.3/.pt,
                position = position_dodge(width = 0.4), size=1/.pt)+
  scale_y_continuous(expand = c(0,0),lim=c(0,max(h2$h2+h2$SE)*1.1))+
  xlab("Genomic regions")+
  theme_Publication3()+
  ylab(expression(paste(h^{2},' for fertility')))


ggsave(filename="Figures/h2.pdf", plot=p, height=1.22, width=1.12, dpi=600)
ggsave(filename="Figures/h2.png", plot=p, height=1.22, width=1.12, dpi=600)
 