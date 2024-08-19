source("./basics_TP.R")
library(gridExtra)
load(file='/Users/tomparee/Documents/Documents - MacBook Pro de tom/Data/rec1EE/outcrossing&outliers/male_frequencies.Rdata')

# common generations 
gcount = unique(subset(male_freq, method=="count")$generation)
gtracker = unique(subset(male_freq, method=="wormtracker")$generation)
common = gcount[which(gcount %in% gtracker)]

data = subset(male_freq, generation %in% common & env == "standard")
data$id = paste0(data$pop, '_', data$generation)

summary(lm(f_male~method + id, data))

# Pearson correlation
data.cast = reshape2::dcast( data , id~method, value.var = "f_male")
cor.test(data.cast$count, data.cast$wormtracker, use = 'complete.obs')




###### RESULTS:
#### In High salt
### Regression:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        0.296580   0.035633   8.323 2.42e-12 ***
# methodwormtracker -0.005921   0.008018  -0.738 0.462517  

### Pearson:.7468621 (P-value)
# p-value = 4.13e-15


#### In NGM
## Regression:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        0.454029   0.023507  19.315  < 2e-16 ***
# methodwormtracker -0.070243   0.007528  -9.330 2.93e-11 ***

## Pearson:
# 0.659137
# p-value = 6.797e-06




##########################
#### Test salt vs env ####
##########################

meanout = aggregate(outcrossing_rate~env+pop, subset(male_freq, method=="wormtracker"), mean,na.rm=T)
sdout.gen = aggregate(outcrossing_rate~env+pop, subset(male_freq, method=="wormtracker"), sd,na.rm=T)
sdout.pop = aggregate(outcrossing_rate~generation+env, subset(male_freq, method=="wormtracker"), sd,na.rm=T)

# Difference in mean outcrossing rate between treatment
t.test(meanout$outcrossing_rate[meanout$env=='standard'],
       meanout$outcrossing_rate[meanout$env=='salt'])

# Difference in mean SD between generation for a pop
t.test(sdout.gen$outcrossing_rate[sdout.gen$env=='standard'], 
       sdout.gen$outcrossing_rate[sdout.gen$env=='salt'])

# Difference in mean SD between pop within a generation
t.test(sdout.pop$outcrossing_rate[sdout.pop$env=='standard'], 
       sdout.pop$outcrossing_rate[sdout.pop$env=='salt'])

#library(car)
#leveneTest(outcrossing_rate ~ env * as.factor(generation), data=subset(male_freq, method=="wormtracker"))

##






##########################
#### Test wt vs mut ####
##########################

meanout = aggregate(outcrossing_rate~rec+pop, subset(male_freq, method=="wormtracker" & env == 'salt'), mean,na.rm=T)


# Difference rec-1 wt and mut
t.test(meanout$outcrossing_rate[meanout$rec=='wt'],
       meanout$outcrossing_rate[meanout$rec=='mut'])

summary(lm(outcrossing_rate~pop+as.factor(generation), data = subset(male_freq, method=="wormtracker" & env == 'salt')))

#library(car)
#leveneTest(outcrossing_rate ~ env * as.factor(generation), data=subset(male_freq, method=="wormtracker"))




library(ggplot2)
library(ggh4x)

mutcolor <- colorRampPalette(c("#F4BB6C", "#9D6329"))
wtcolor<- colorRampPalette(c("lightblue2", "#1F476E"))

agg.male_freq = aggregate(rec~pop+env, data=male_freq, function(x){x[1]})
agg.male_freq$color = NA
agg.male_freq$color[agg.male_freq$env=="salt" & agg.male_freq$rec=="wt"] = wtcolor(8)
agg.male_freq$color[agg.male_freq$env=="salt" & agg.male_freq$rec=="mut"] = mutcolor(8)
agg.male_freq$color[agg.male_freq$env=="standard" & agg.male_freq$rec=="wt"] = wtcolor(4)
agg.male_freq$color[agg.male_freq$env=="standard" & agg.male_freq$rec=="mut"] = mutcolor(4)

p1=ggplot(male_freq, aes(generation, f_male, linetype=method, color=rec))+
  theme_Publication3()+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  geom_point(shape=1, size= 0.8/.pt)+
  geom_line(size= 1.3/.pt)+
  facet_wrap2(~pop, ncol=8, scales='free')+
  scale_color_manual(values = c(orangeMut, blueWT))+
  theme(legend.position = "none")+
  xlab("Generation")+
  ylab("Male frequency")+
  ylim(range(male_freq$f_male))+
  xlim(range(male_freq$generation))+
  geom_text(data=data.frame(pop = c("SHR4", "SLR3"), method = "count", rec='wt'), 
            aes(x=22, y=0.65, label = "*"), 
            color="black", size= 11/.pt)


p2=ggplot(subset(male_freq, method=="wormtracker"), aes(generation-1,outcrossing_rate, color=pop))+
  geom_line(aes(group=pop), linetype='dashed', alpha=0.7, size=1/.pt)+
  facet_wrap(~factor(env, levels = c("standard", "salt"), labels = c("Domestication", "Novel")), nrow=2, scales='free')+
  theme_Publication3()+
  theme(strip.text = element_text(size=7))+
  scale_color_manual(values = c(agg.male_freq$color, orangeMut, blueWT), breaks = c(agg.male_freq$pop, "mut", "wt"))+
  ylim(0,1)+stat_summary(geom='line', aes(color=rec), size=1/.pt)+
  stat_summary(aes(color=rec), size=0.3/.pt)+ theme(legend.position="none")+
  xlab("Generation")+ylab("Outcrossing rate")


p=grid.arrange(p1,p2+theme(plot.margin = unit(c(10,10,10,30), 'pt') ), ncol=2, widths = c(2.5,1))

ggsave(filename="/Users/tomparee/Desktop/outcrossing.png", plot=p, height=6.5*weirdFigureRatio, width=18.3*weirdFigureRatio, units = "cm", dpi=1000)
ggsave(filename="/Users/tomparee/Desktop/outcrossing.pdf", plot=p, height=6.5*weirdFigureRatio, width=18.3*weirdFigureRatio, units = "cm", dpi=1000)

