####################################
#### MANATHAN PLOT #################
####################################
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")
library(ggplot2)

load(file = "candidateSelection/candidate_snps.RData")
load("./genotype/Rhom_poolseq.RData")


permutated_min_pval <- read.csv("candidateSelection/permutated_min_pval.txt", row.names=NULL, sep="")
permutated_min_pval = permutated_min_pval[,-1]
#write.table(permutated_min_pval,file="./Data/Experimental_Evolution_R/EE_Rhom/poolseq/candidate_loci/permutated_min_pval_1to600.txt")
#neutral <- read.csv("./Data/Experimental_Evolution_R/EE_Rhom/poolseq/candidate_loci/permutated_min_pval.txt", sep="")
neutral = permutated_min_pval
threshold = apply(neutral, 2, function(x){
  x = x[order(x)]
  threshold = 0.05
  x[length(x)*threshold]
})

threshold = data.frame(threshold = threshold , test = names(threshold ))

marker = candidate_snps$marker


wtest = c('pstandard', 'psalt',
          'pstandardxrec','psaltxrec')

candidate_snps = candidate_snps[,wtest]

significant = apply(candidate_snps, 1, function(x){sum(x<threshold$threshold)})

candidate_snps = reshape2::melt(candidate_snps)
names(candidate_snps) = c('test', 'pval')
candidate_snps$marker = rep(marker, length(wtest))
candidate_snps$chrom = data.table::tstrsplit(candidate_snps$marker, "_")[[1]]
candidate_snps$pos = as.numeric(data.table::tstrsplit(candidate_snps$marker, "_")[[2]])



candidate_snps= do.call(rbind, lapply(split(candidate_snps, candidate_snps$test), function(x){
  
  test = x$test[1]
  threshold = threshold$threshold[threshold$test==test]
  
  x$significant = ifelse(x$pval<threshold, "yes", "no")
  
  x
  
}))

candidate_snps = subset(candidate_snps, test %in% c('pstandard', 'psalt',
                                        'pstandardxrec', 'psaltxrec'))
candidate_snps$env = tstrsplit(candidate_snps$test, 'x',fixed = TRUE)[[1]]
candidate_snps$logp = -log10(candidate_snps$pval) + 0.5
candidate_snps$logp[!grepl('x',candidate_snps$test)] = -candidate_snps$logp[!grepl('x',candidate_snps$test)]




df.threshold = threshold
df.threshold$env = ifelse(grepl("standard",df.threshold$test), "pstandard", "psalt")
df.threshold$logth = -log10(df.threshold$threshold)

df.threshold$logth = ifelse(grepl("xrec",df.threshold$test), 1, -1)*df.threshold$logth 


candidate_snps = merge(candidate_snps, data.frame(chrom=unique(candidate_snps$chrom), odd = ifelse((1:length(unique(candidate_snps$chrom))) %% 2 == 0, 0, 1)))
candidate_snps$pos = as.numeric(data.table::tstrsplit(candidate_snps$marker, "_")[[2]])
maxpos = aggregate(pos~chrom, candidate_snps, max)
candidate_snps = merge(candidate_snps, data.frame(chrom=maxpos$chrom,pos_toadd=c(0, cumsum(maxpos[1:5,]$pos))))
candidate_snps$pos2 = candidate_snps$pos + candidate_snps$pos_toadd
candidate_snps$test_type = ifelse(grepl("xrec",candidate_snps$test), "genxrec", "gen")
chromticks = aggregate(pos2~chrom, candidate_snps, mean)

#"#43B547","#288E2B"
#"#5B86BD","#73A8EE"
#"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
# p=ggplot()+
#   facet_wrap(~factor(env, level=c('pstandard', 'psalt'), labels = c("Domestication", "High salt")), nrow=2, scale = "free")+
#   theme_Publication3()+theme(panel.spacing = unit(0, "pt"), strip.text.x = element_text(size = 7))+
#   theme(axis.text.y = element_text(size = 6),axis.text.x = element_text(size = 6),axis.title.y = element_text(size = 7), axis.title.x = element_blank(), legend.justification  = "top")+
#   geom_point(data = candidate_snps, aes(pos2/1e6, logp, color = paste0(test_type, odd)), size=0.3/.pt)+
#   scale_color_manual(values = c(genxrec1 = "#3C5488FF",genxrec0="#8491B4FF", gen1="#00A087FF",gen0="#91D1C2FF"),
#                      breaks = c("genxrec1","gen1"), name = "", labels = c("~ generation x rec-1", "~ generation"))+
#   scale_y_continuous(breaks = c(-10.5, -5.5,0,5.5), labels = c(10,5,0, 5), expand = c(0.05,0))+
#   geom_hline(data = df.threshold, aes(yintercept=logth), size=1/.pt, color = orangeMut)+
#   scale_x_continuous(expand = c(0,0), breaks = chromticks$pos2/1e6, labels = chromticks$chrom)+
#   ylab("-log10(P-value)")+
#   guides(colour = guide_legend(override.aes = list(size=4)))+
#   theme(legend.position = "bottom", legend.key.size = unit(1, 'pt'))
# 

p=ggplot()+
  facet_wrap(~factor(env, level=c('pstandard', 'psalt'), labels = c("Domestication", "High salt")), nrow=1, scales = 'free')+
  theme_Publication3()+theme(panel.spacing = unit(0, "pt"), strip.text.x = element_text(size = 7))+
  theme(axis.text.y = element_text(size = 6),axis.text.x = element_text(size = 6),axis.title.y = element_text(size = 7), axis.title.x = element_blank(), legend.justification  = "top")+
  geom_point(data = candidate_snps, aes(pos2/1e6, logp, color = paste0(test_type, odd)), size=0.1/.pt)+
  scale_color_manual(values = c(genxrec1 = "#3C5488FF",genxrec0="#8491B4FF", gen1="#00A087FF",gen0="#91D1C2FF"),
                     breaks = c("genxrec1","gen1"), name = "", labels = c("~ generation x rec-1", "~ generation"))+
  scale_y_continuous(breaks = seq(-20,20,5), labels =  abs(seq(-20,20,5)), expand = c(0.05,0), lim=range(candidate_snps$logp))+
  geom_hline(data = df.threshold, aes(yintercept=logth), size=1/.pt, color = orangeMut)+
  scale_x_continuous(expand = c(0,0), breaks = chromticks$pos2/1e6, labels = chromticks$chrom)+
  ylab("-log10(P-value)")+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.position = "bottom", legend.key.size = unit(1, 'pt'))


ggsave(filename="Figures/manhattan.png", plot=p, height=1.5, width=3.3, dpi=600)
ggsave(filename="/Figures/manhattan.pdf", plot=p,height=1.5, width=3.3, dpi=600)



# p=ggplot()+
#   facet_wrap(~factor(env, level=c('pstandard', 'psalt'), labels = c("Domestication", "Novel")), nrow=1)+
#   theme_Publication3()+theme(panel.spacing = unit(15, "pt"), strip.text.x = element_text(size = 7))+
#   theme(axis.text.y = element_text(size = 6),axis.text.x = element_text(size = 6),axis.title.y = element_text(size = 7), axis.title.x = element_blank(), legend.justification  = "top")+
#   geom_point(data = candidate_snps, aes(pos2/1e6, logp, color = paste0(test_type, odd)), size=0.25/.pt)+
#   scale_color_manual(values = c(genxrec1 = "#3C5488FF",genxrec0="#8491B4FF", gen1="#00A087FF",gen0="#91D1C2FF"),
#                      breaks = c("genxrec1","gen1"), name = "", labels = c("~ generation x rec-1", "~ generation"))+
#   scale_y_continuous(breaks = c(-10.5, -5.5,0,5.5), labels = c(10,5,0, 5), expand = c(0.05,0))+
#   geom_hline(data = df.threshold, aes(yintercept=logth), size=1/.pt, color = orangeMut)+
#   scale_x_continuous(expand = c(0,0), breaks = chromticks$pos2/1e6, labels = chromticks$chrom)+
#   ylab("-log10(P-value)")+
#   guides(colour = guide_legend(override.aes = list(size=4)))+
#   theme(legend.position = "bottom", legend.key.size = unit(1, 'pt'))
# 
# 
# ggsave(filename="/Users/tomparee/Desktop/Figures/manhattan_h.png", plot=p, height=2.5, width=4, dpi=600)
# ggsave(filename="/Users/tomparee/Desktop//Figures/manhattan_h.pdf", plot=p, height=2.5, width=4, dpi=600)
# 
# 
# 
# 
# 
# 


########################################################################
######## MANATHAN GEN X ENV


load("./candidateSelection/candidate_snps.RData")
load("./genotype/Rhom_poolseq.RData")

permutated_min_pval <- read.csv("candidateSelection/permutated_min_pval_sGenxenv.txt", row.names=NULL, sep="")
permutated_min_pval = permutated_min_pval[,-1]

#write.table(permutated_min_pval,file="./Data/Experimental_Evolution_R/EE_Rhom/poolseq/candidate_loci/permutated_min_pval_1to600.txt")
#neutral <- read.csv("./Data/Experimental_Evolution_R/EE_Rhom/poolseq/candidate_loci/permutated_min_pval.txt", sep="")
neutral = permutated_min_pval
threshold = apply(neutral, 2, function(x){
  x = x[order(x)]
  threshold = 0.05
  x[length(x)*threshold]
})

threshold = data.frame(threshold = threshold , test = names(threshold ))

marker = candidate_snps$marker

wtest=c("sGenxenv")

candidate_snps = candidate_snps[,wtest]
candidate_snps = as.data.frame(matrix(candidate_snps, ncol=1))

#ignificant = candidate_snps<threshold$threshold[threshold$test = sGenxenv]


candidate_snps$marker = rep(marker, length(wtest))
candidate_snps$chrom = data.table::tstrsplit(candidate_snps$marker, "_")[[1]]
candidate_snps$pos = as.numeric(data.table::tstrsplit(candidate_snps$marker, "_")[[2]])
candidate_snps$logp = -log(candidate_snps$V1, base=10)
#candidate_snps$significant = candidate_snps[,1]<threshold$threshold


maxpos = aggregate(pos~chrom, candidate_snps, max)
candidate_snps = merge(candidate_snps, data.frame(chrom=maxpos$chrom,pos_toadd=c(0, cumsum(maxpos[1:5,]$pos))))
candidate_snps$pos2 = candidate_snps$pos + candidate_snps$pos_toadd
chromticks = aggregate(pos2~chrom, candidate_snps, function(x){(min(x)+max(x))/2})
candidate_snps$odd = ifelse(as.numeric(romtonum(candidate_snps$chrom)) %% 2 == 0, 0, 1)


data.frame(chrom=maxpos$chrom,pos_toadd=c(0, cumsum(maxpos[1:5,]$pos)))

df.threshold = threshold
df.threshold$env = ifelse(grepl("standard",df.threshold$test), "pstandard", "psalt")
df.threshold$logth = -log10(df.threshold$threshold)
df.threshold = subset(df.threshold, test =="sGenxenv")



p=ggplot()+
  theme_Publication3()+
  geom_point(data = candidate_snps, aes(pos2/1e6, logp, color=as.factor(odd)), size=0.4/.pt)+
  geom_hline(data = df.threshold, aes(yintercept=logth), size=1/.pt, color = orangeMut)+
  ylab("-log10(P-value)")+xlab("")+
  scale_color_manual(values = c("black", "darkgrey"))+
  scale_x_continuous(expand = c(0,0), breaks = chromticks$pos2/1e6, labels = chromticks$chrom)+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = "none")


ggsave(filename="/Users/tomparee/Desktop/manhattan_GxENV.png", plot=p, height=1.5, width=6, dpi=600)
ggsave(filename="/Users/tomparee/Desktop/manhattan_GxENV.pdf", plot=p, height=1.5, width=6, dpi=600)






