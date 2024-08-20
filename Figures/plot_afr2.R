
library(data.table)
library(ggplot2)

PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

load(file = "observed_r2/afr2.Rdata")


#afr2.mean = do.call(rbind, lapply(split( afr2, paste0( afr2$rec, afr2$rrchange)), function(x){
#  
#  
#  #x = split( afr2, paste0( afr2$rec, afr2$rrchange))[[1]]
#  
#  out = apply(x[,grepl("mean.r2|boot",colnames(x))],2, function(y){weighted.mean(y, x[,"nobs"])})
#  
#  out = as.data.frame(matrix(out, nrow=1))
#  out = cbind(rrchange = x[1,"rrchange"],x[1,"rec"],sum(x$nobs),out)
#  colnames(out)=c("rrchange.region", "rec","nobs", colnames(x)[grepl("mean.r2|boot",colnames(x))])
#  out
#}))

afr2.mean = do.call(rbind, lapply(split( afr2, paste0( afr2$rec, afr2$rrchange)), function(x){
  
  out = apply(x[,grepl("mean.r2|boot",colnames(x))],2, function(y){weighted.mean(y, x[,8])})
  out = setNames(as.data.frame(matrix(out, ncol=1)), "mean.r2")
  out$type = c("observed", rep("bootstrap",nrow(out)-1))
  out$n = sum(x$nobs)
  out$rrchange.region= x[1,"rrchange"]
  out$rec = x[1,"rec"]
  out
}))



bootsub = subset(r2.mean, type == "bootstrap")

CI= do.call(rbind, lapply(split(bootsub, paste0(bootsub$rrchange.region,bootsub$rec, bootsub$n)), function(x){
  qx = setNames(as.data.frame(matrix(quantile(x$mean.r2, c(0.005, 0.025, 0.975, 0.995)), nrow=1)), c("low.CI.0.99", "low.CI.0.95", "high.CI.0.95","high.CI.0.99"))
  cbind(x[1,c("n", "rrchange.region", "rec")], qx)
}))


afr2.mean = merge(subset(afr2.mean, type == "observed"), CI)
afr2.mean = afr2.mean[,-which(colnames(afr2.mean)=="type")]








afr2.mean$rrchange.region = factor(afr2.mean$rrchange.region, levels = c("genome-wide","1", "-1"))

ph=ggplot(afr2.mean)+theme_Publication3()+
  geom_errorbarh(aes(y=as.factor(rrchange.region), xmin = low.CI.0.95, xmax = high.CI.0.95, color =rec), height = 0.7/.pt, size= 1/.pt)+
  geom_point(aes(y=as.factor(rrchange.region), x=mean.r2, color = rec), size =1.2/.pt)+
  geom_errorbarh(data=subset(afr2.mean, rrchange.region=="genome-wide"), aes(y=as.factor(rrchange.region), xmin = low.CI.0.95, xmax = high.CI.0.95, color =rec), height = 0.8/.pt, size= 1.4/.pt)+
  geom_point(data=subset(afr2.mean, rrchange.region=="genome-wide"), aes(y=as.factor(rrchange.region), x=mean.r2, color = rec), size =1.7/.pt)+
  scale_y_discrete(breaks = c("genome-wide","1", "-1"), labels = c("Genome-wide", "Increased rr","Decreased rr"))+
  scale_color_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut), name="", labels = c("Wild-type", "Mutant"))+
  ylab("")+xlab("Allele frequencies correlation")


afr2.mean$rrchange.region = factor(afr2.mean$rrchange.region, levels = c("1", "-1","genome-wide"))
pv=ggplot(afr2.mean)+theme_Publication3()+
  geom_errorbar(aes(x=as.factor(rrchange.region), ymin = low.CI.0.95, ymax = high.CI.0.95, color =rec), width = 0.7/.pt, size= 1/.pt)+
  geom_point(aes(x=as.factor(rrchange.region), y=mean.r2, color = rec), size =1.2/.pt)+
  geom_errorbar(data=subset(afr2.mean, rrchange.region=="genome-wide"), aes(x=as.factor(rrchange.region), ymin = low.CI.0.95, ymax = high.CI.0.95, color =rec), width = 0.8/.pt, size= 1.4/.pt)+
  geom_point(data=subset(afr2.mean, rrchange.region=="genome-wide"), aes(x=as.factor(rrchange.region), y=mean.r2, color = rec), size =1.7/.pt)+
  scale_x_discrete(breaks = c("genome-wide","1", "-1"), labels = c("Genome-wide", "Increased rr","Decreased rr"))+
  scale_color_manual(breaks = c("wt", "mut"), values = c(blueWT, orangeMut), name="", labels = c("Wild-type", "Mutant"))+
  ylab("Allele frequencies correlatio")+xlab("")




ggsave(filename="Figures/afr2.png", plot=ph, height=1, width=2.6, dpi=600)
ggsave(filename="Figures/afr2.pdf", plot=ph, height=1, width=2.6, dpi=600)





