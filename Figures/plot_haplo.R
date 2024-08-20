PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")
source("haplotype/function_findhaplo_in_RILS.R")
snps = as.data.frame(fread("genotype/chrI_SMR_inbredlines_snps_CeMEEv2_ws245.csv.gz"))
load(file="haplotype/RILs_haplotypes_foundersOnly.Rdata")



load(file="haplotype/RILs_haplotypes_founders&EEV1403_resolved.Rdata")

showEEV1403 = F

if(!showEEV1403){
  haplotype1 = haplo
  load(file="haplotype/RILs_haplotypes_foundersOnly_resolved.Rdata")
  haploEEV1403 = subset(haplo, rilname == "EEV1401")
  
  
  haplotype = do.call(rbind, lapply(1:nrow(haplotype1), function(i){
    if(i %% 1000 == 0){print(i)}
    x = haplotype1[i,]
    if(x[,"EEV1403"]==1){
      whichsnp1 = x[,"whichsnp1"]
      whichsnp2 = x[,"whichsnp2"]
      
      x2 = haploEEV1403[haploEEV1403[,"whichsnp1"] <= whichsnp2 & haploEEV1403[,"whichsnp2"] >= whichsnp1,]
      x2[1,"whichsnp1"] = whichsnp1
      x2[nrow(x2),"whichsnp2"] = whichsnp2
      x2[,"rilname"] = x[,"rilname"]
    }else{
      x2 = x[,-which(colnames(x)=="EEV1403")]
    }
    
    return(x2)
    
  }))
  
  haplo = rbind(haplotype,haploEEV1403)
  haplo=haplo[!is.na(haplo$AB1),]
}







orderedRILS = sort(unique(haplo$rilname))
orderedRILS=c(orderedRILS[grepl("mut",orderedRILS)], orderedRILS[grepl("wt",orderedRILS)], orderedRILS[grepl("EEV1401",orderedRILS)])
foundernames = colnames(haplo)[3:18]
haplop = haploPlotFromat(Rilshaplotypes=haplo,
                        widthRILs = 0.85,
                        foundercolnames=foundernames,
                        founderscolors = NULL,na.color=NULL, info=NULL,
                        orderedRils = orderedRILS)

haplop$pos1 = snps$pos[haplop$whichsnp1]
haplop$pos2 = snps$pos[haplop$whichsnp2]

haplop$y1[grepl(haplo)]


colortable = data.frame(founders = foundernames,
                        haplocolor = c("black","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF",
                                             "#E18726", "#332F92", "#FFDC91", "#008886", "#90C98C"))
ggplot(haplop)+theme_Publication3()+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6, ymin=y1, ymax=y2,fill=founder), color = NA)+
  scale_fill_manual(values = colortable$haplocolor, breaks = colortable$founder,labels = colortable$founder)+
  scale_x_continuous(expand = c(0,0), breaks = seq(1,16,1), lim = c(-0.2, 15.2))+
  #scale_y_continuous(expand = c(0,0), breaks = c(27,77)+7, labels = c("Mutant","Wild-type"),lim = c(0, max(haplop$y)+9))+
  #theme(axis.text.y = element_text(angle = 90, size = 7),
  #      axis.line.y = element_line(color = "white"),
  #      axis.ticks.y = element_line(color = "white"),
  #      axis.title.y = element_blank(),
  #      axis.text.x = element_text(size = 7),
  #      legend.key.size = unit(10, 'pt'))+
  ylab("")+
  #geom_segment(x=-0.1,xend=-0.1,y=0.5,yend=50.5, linewidth = 1/.pt)+
  #geom_segment(x=-0.1,xend=-0.1,y=52.5,yend=102, linewidth = 1/.pt)+
  #geom_point(aes(x= 717819/1e6, y= max(haplop$y)+2), shape = 25, fill = "black", size = 4.5/.pt)+
  #geom_text(aes(x= 717819/1e6, y= max(haplop$y)+6, label = 'rec-1'), size = 7/.pt, color = "black", fontface = "italic")+#fontface=3
  xlab("Physical size (Mb)")

