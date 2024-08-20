##################################################################################
##################### plot SNV density & delta R  #####################################
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggsci)
source("./rec1EE_utils.R")

load("./recombinationMaps/recmaps.Rdata")
load("./recombinationMaps/deltaR.Rdata")
load("genotype/Rhom_poolseq.RData")
##########################################################
##########################################################


tarchrom = c("I", "II", "III", "IV", "V", "X") # chromosomes we want to study, here all

#############################################@
# Calculate SNV density along the chromosome

densityWinsize = 1e5
density = do.call(rbind, lapply(split(snps, snps$chrom), function(x){
  
  windows = get.win(size =max(x$pos), winsize = densityWinsize, minsize = densityWinsize/2)
  
  do.call(rbind, lapply(windows, function(win){
    pos1 = min(win)
    pos2 = max(win)
    xwin = subset(x, pos>=pos1 & pos<pos2)
    density = nrow(xwin)/(pos2-pos1+1)
    data.frame(chrom = x$chrom[1],pos1 = pos1, pos2=pos2, snpdensity = density)
  }))
  
}))

##############
###PLOT:  ####

##########################################################
### SNV density and delta R (FIG 2 AB): color version ####

plotDeltaR=ggplot(recchange)+theme_Publication3()+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0, 0),breaks = seq(0,25,4))+
  scale_y_continuous(expand = c(0, 0),breaks = c(-8,-4,0,4,8))+
  geom_point(data=data.frame(x=0,y=0, chrom = romtonum(tarchrom)),aes(x,y), color=NA)+
  geom_hline(yintercept = 0, size = 1/.pt, linetype = 'dashed')+
  geom_step(aes(pos1/1e6, rrdiff), size=1/.pt,color = "grey55")+
  geom_segment(aes(x=pos1/1e6, xend = pos2/1e6, y=rrdiff, yend = rrdiff), size=2/.pt,color = "grey55")+
  geom_point(data=data.frame(chrom=1, pos = recpos/1e6), aes(pos, max(recchange$rrdiff)*1.2), shape=25, size=3/.pt, fill='black')+
  geom_text(data=data.frame(chrom=1, pos = recpos/1e6), aes(pos, max(recchange$rrdiff)*1.7, label = "rec-1"), size = 7/.pt, color='black',fontface="italic")+
  geom_line(data=snps, aes(pos/1e6, rrdiff, color=chrom), size=0.5)+
  ylab("\u{0394}r" )+
  xlab("Physical position (Mb)")+
  coord_cartesian(clip="off", ylim = c(min(recchange$rrdiff)*1.1,max(recchange$rrdiff)))+
  facet_grid(.~numtorom(chrom), scales = "free_x", space = "free")+
  theme(panel.spacing = unit(10, "pt"))+
  theme(axis.text = element_text(size=6.5))+
  scale_color_npg()+theme(legend.position="none")

plotSNPdensity=ggplot(subset(density, chrom %in% tarchrom))+theme_Publication3()+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6,ymin= 0, ymax = snpdensity*1e3, fill=chrom), alpha=0.8)+
  geom_step(aes(x=pos1/1e6,y = snpdensity*1e3), size = 1/.pt)+
  scale_y_continuous(expand = c(0, 0))+
  geom_point(data=data.frame(x=0,y=0), aes(x,y), color=NA)+
  #scale_x_continuous(expand = c(0,0),breaks = seq(0,25,2), lim = c(0,max(subset(density, chrom=='I')$pos2)/1e6))+
  scale_x_continuous(expand = c(0,0),breaks = seq(0,25,4))+
  ylab("SNV density")+
  xlab("Physical position (Mb)")+
  facet_grid(.~chrom, scales = "free_x", space = "free")+
  theme(panel.spacing = unit(10, "pt"))+
  theme(axis.text = element_text(size=6.5))+
  scale_fill_npg()+theme(legend.position="none")




p=grid.arrange(plotDeltaR+theme(axis.title.x = element_blank(),plot.margin = unit(c(30,10,5,11), "pt")),
               plotSNPdensity+theme(plot.margin = unit(c(5,10,10,10), "pt"),strip.background = element_blank(), strip.text = element_blank()),
               heights = c(1,0.8))

ggsave(filename="Figures/rrdiff_color.png", plot=p, height=2.48, width=5.4, dpi=600)
ggsave(filename="Figures/rrdiff_color.pdf", plot=p,height=2.48, width=5.4, dpi=600)


#######################################################################
### SNV density and delta R: grey & black version #####################
### + recombination maps (for HT) #####################################


plotMarey=ggplot(recmaps)+theme_Publication3()+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0, 0),breaks = seq(0,25,4))+
  #scale_y_continuous(expand = c(0, 0))+
  geom_line(aes(pos/1e6, genetic, color=rec), size=0.5)+
  geom_text(data=data.frame(chrom=1, pos = recpos/1e6), aes(pos, 58, label = "rec-1"), size = 7/.pt, color='black',fontface="italic")+
  geom_point(data=data.frame(chrom=1, pos = recpos/1e6), aes(pos, 51), shape=25, size=3/.pt, fill='black')+
  scale_color_manual(values = c(orangeMut,blueWT))+
  ylab("Genetic distance (cM)" )+
  xlab("Physical position (Mb)")+
  coord_cartesian(clip="off",ylim=c(0,50),expand=0)+
  facet_grid(.~numtorom(chrom), scales = "free_x", space = "free")+
  theme(panel.spacing = unit(10, "pt"))+
  theme(axis.text = element_text(size=6.5))+
  theme(legend.position="none")

plotDeltaR2=ggplot(recchange)+theme_Publication3()+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0, 0),breaks = seq(0,25,4))+
  scale_y_continuous(expand = c(0, 0),breaks = c(-8,-4,0,4,8))+
  geom_point(data=data.frame(x=0,y=0, chrom = romtonum(tarchrom)),aes(x,y), color=NA)+
  geom_hline(yintercept = 0, size = 1/.pt, linetype = 'dashed')+
  geom_step(aes(pos1/1e6, rrdiff), size=1/.pt,color = "grey55")+
  geom_segment(aes(x=pos1/1e6, xend = pos2/1e6, y=rrdiff, yend = rrdiff), size=2/.pt,color = "grey55")+
  #geom_point(data=data.frame(chrom=1, pos = recpos/1e6), aes(pos, max(recchange$rrdiff)*1.2), shape=25, size=3/.pt, fill='black')+
  #geom_text(data=data.frame(chrom=1, pos = recpos/1e6), aes(pos, max(recchange$rrdiff)*1.7, label = "rec-1"), size = 7/.pt, color='black')+
  geom_line(data=snps, aes(pos/1e6, rrdiff), size=0.5)+
  ylab("\u{0394}r" )+
  xlab("Physical position (Mb)")+
  #coord_cartesian(clip="off", ylim = c(min(recchange$rrdiff)*1.1,max(recchange$rrdiff)))+
  facet_grid(.~numtorom(chrom), scales = "free_x", space = "free")+
  theme(panel.spacing = unit(10, "pt"))+
  theme(axis.text = element_text(size=6.5))+
  theme(legend.position="none")

plotSNPdensity2=ggplot(subset(density, chrom %in% tarchrom))+theme_Publication3()+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6,ymin= 0, ymax = snpdensity*1e3), alpha=0.8, fill='grey')+
  geom_step(aes(x=pos1/1e6,y = snpdensity*1e3), size = 1/.pt)+
  scale_y_continuous(expand = c(0, 0))+
  geom_point(data=data.frame(x=0,y=0), aes(x,y), color=NA)+
  #scale_x_continuous(expand = c(0,0),breaks = seq(0,25,2), lim = c(0,max(subset(density, chrom=='I')$pos2)/1e6))+
  scale_x_continuous(expand = c(0,0),breaks = seq(0,25,4))+
  ylab("SNV density")+
  xlab("Physical position (Mb)")+
  facet_grid(.~chrom, scales = "free_x", space = "free")+
  theme(panel.spacing = unit(10, "pt"))+
  theme(axis.text = element_text(size=6.5))+
  theme(legend.position="none")

p2=grid.arrange(plotMarey+theme(axis.title.x = element_blank(),plot.margin = unit(c(30,10,5,11), "pt")),
                plotDeltaR2+theme(axis.title.x = element_blank(), plot.margin = unit(c(5,10,10,10), "pt")),
               plotSNPdensity2+theme(plot.margin = unit(c(5,10,10,10), "pt"),strip.background = element_blank(), strip.text = element_blank()),
               heights = c(1,0.8,0.8))


ggsave(filename="Figures/rrdiff_HT.png", plot=p2, height=3.4, width=5.4, dpi=600)
ggsave(filename="Figures/rrdiff_HT.pdf", plot=p2,height=3.4, width=5.4, dpi=600)




##################################################################
### SNV density ~ delta R ########################################

recmaps$chrom = numtorom(recmaps$chrom)
winsize=1e6
winstat = do.call(rbind, lapply(split(recmaps, recmaps$chrom), function(x){
  
  #x = split(recmaps, recmaps$chrom)[[2]]
  chr = x$chrom[1]
  ppos = seq(1,max(x$pos),winsize)
  ppos[length(ppos)] = max(x$pos)
  
  gwt = subset(recmaps, rec=='wt' & chrom == chr)
  gmut = subset(recmaps, rec=='mut' & chrom == chr)
  
  gwt = approx(x = gwt$pos, y = gwt$genetic, xout = ppos)$y
  gmut = approx(x = gmut$pos, y = gmut$genetic, xout = ppos)$y
  dpos = diff(ppos)
  rrwt = diff(gwt)/dpos
  rrmut = diff(gmut)/dpos
  
  density = sapply(1:(length(ppos)-1), function(i){nrow(subset(snps, chrom==chr & pos>= ppos[i] & pos < ppos[i+1]))})/dpos
  
  out = data.frame(chrom=chr,
                   pos1 = ppos[-length(ppos)],
                   pos2 = ppos[-1],
                   rrwt = rrwt,
                   rrmut = rrmut,
                   rrdiff = rrmut-rrwt,
                   density = density,
                   nsnp = density*dpos)
  
  out=out[order(out$rrdiff),]
  out$cumdensity = cumsum(out$density)
  out
  
}))





p3=ggplot(winstat, aes(rrdiff*1e6,density*1000))+theme_Publication3()+
  geom_point(size =  1.5/.pt, aes(color = chrom), alpha =0.5)+
  geom_smooth(method = "lm", se = F, aes(group = chrom, color = chrom), size =  1.7/.pt)+
  geom_smooth(method = "lm", se = F, size = 1/.pt, linetype = "dashed", color = "black",fullrange=TRUE,size =  2/.pt)+
  scale_x_continuous(expand = c(0, 0), limits = max(abs(winstat$rrdiff*1e6))*c(-1.1,1.1))+
  scale_color_npg(name = "")+
  xlab("\u{0394}r" )+
  ylab("SNV density")+
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6))+
  guides(fill=guide_legend(keywidth=10,keyheight=10,default.unit="inch"))



ggsave(filename="Figures/density~rr.pdf", plot=p3, height=1.22, width=1.77, dpi=600)
ggsave(filename="Figures/density~rr.png", plot=p3, height=1.22, width=1.77, dpi=600)

