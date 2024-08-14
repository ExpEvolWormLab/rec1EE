
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")
library(emmeans)
library(data.table)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)

fitness = as.data.frame(fread("adaptation_fitness/fitness.txt"))



####################################
### ANCESTRAL FITNESS ###############

ancestor = subset(fitness, generation == "G0")
ancestor$rec = factor(ancestor$rec, levels = c("wt", "mut"))
ancestor = do.call(rbind, lapply(split(ancestor, ancestor$env), function(x){
  e = x$env[1]
  wtfit = subset(mean.ancestral.fitness, env==e & rec=="wt")$rel.fitness
  x$rel.fitness =  x$rel.fitness/wtfit
  x
}))

mean.ancestral.fitness = do.call(rbind, lapply(split(fitness, fitness$env), function(x){
  model <- lm(rel.fitness~block+rec, data = subset(x,generation=="G0"))
  value = emmeans(model, specs = pairwise ~ rec, level= 0.95)
  value = as.data.frame(value$emmeans)
  colnames(value)[colnames(value)=="emmean"]="rel.fitness"
  value$env =x$env[1]
  value
}))


mean.ancestral.fitness2 = do.call(rbind,lapply(split(mean.ancestral.fitness, mean.ancestral.fitness$env), function(x){
  wtfit = subset(x, rec=='wt')$rel.fitness
  x$rel.fitness =  x$rel.fitness/wtfit
  x$lower.CL = x$lower.CL/wtfit
  x$upper.CL = x$upper.CL/wtfit
  x
}))


mean.ancestral.fitness$rec = factor(mean.ancestral.fitness$rec, levels = c('wt','mut'))
mean.ancestral.fitness$env = factor(mean.ancestral.fitness$env, levels = c('ngm','salt'))





#################################
### FITNESS GAIN ###############

# Model average wt and mut per env
mean.fitness.gain.rec = do.call(rbind, lapply(split(fitness, fitness$env), function(x){
  model <- lmer(fitness.gain.rel~(1|pop)+(1|block)+rec, data = subset(x, generation=="G40"))
  
  #fitness.gain.wt = summary(model)$coef[1,1]
  #fitness.gain.mut = fitness.gain.wt + summary(model)$coef[2,1]
  
  value = emmeans(model, specs = pairwise ~ rec, type = "response", level=0.95)
  value = as.data.frame(value$emmeans)
  colnames(value)[colnames(value)=="emmean"]="fitness.gain.rel"
  value$env = x$env[1]
  value
}))


# Model average per pop
mean.fitness.gain.pop = do.call(rbind, lapply(split(subset(fitness, generation=="G40"), subset(fitness, generation=="G40")$pop), function(x){
  #x=split(subset(fitness, generation=="G40"), subset(fitness, generation=="G40")$pop)[[1]]
  model <- lmer(fitness.gain.rel~(1|block), data = x)
  mx = summary(model)$coef[1,1]
  SEx =  summary(model)$coef[1,2]
  value = data.frame(pop=x$pop[1], fitness.gain.rel=mx, SE=SEx, lower.CL = mx-(1.96*SEx), upper.CL = mx+(1.96*SEx))
  value
}))

# model <- lmer(fitness.gain.rel~(1|block)+pop, data = subset(fitness, generation=="G40"))
# 
# fitness.gain.wt = summary(model)$coef[1,1]
# fitness.gain.mut = fitness.gain.wt + summary(model)$coef[2,1]
# 
# value = emmeans(model, specs = pairwise ~ pop, type = "response", level=0.95)
# value = as.data.frame(value$emmeans)
# colnames(value)[colnames(value)=="emmean"]="fitness.gain.rel"
# 

mean.fitness.gain.pop = assign.infopop.Rhom(mean.fitness.gain.pop)

mean.fitness.gain.pop = mean.fitness.gain.pop[order(mean.fitness.gain.pop$env, mean.fitness.gain.pop$rec,mean.fitness.gain.pop$fitness.gain.rel, decreasing = c(F,F,F)),]
mean.fitness.gain.pop$pop_ordered = as.numeric(factor(mean.fitness.gain.pop$pop, levels = mean.fitness.gain.pop$pop, labels = 1:length(unique(mean.fitness.gain.pop$pop))))
mean.fitness.gain.pop$pop_ordered[mean.fitness.gain.pop$rec=="mut"] = mean.fitness.gain.pop$pop_ordered[mean.fitness.gain.pop$rec=="mut"]+4
xticks = sort(aggregate(pop_ordered~rec+env,unique(mean.fitness.gain.pop[,c("env", "rec", "pop_ordered")]), mean)$pop_ordered)

mean.fitness.gain.pop = do.call(rbind, lapply(split(mean.fitness.gain.pop, paste0(mean.fitness.gain.pop$env, mean.fitness.gain.pop$rec)), function(x){
  
  #x=split(adapted, paste0(adapted$env, adapted$rec))[[1]]
  
  po = unique(x$pop_ordered)
  pr = sample(po, length(po))
  x$pop_random = pr[match(x$pop_ordered, po)]
  x
  
}))


adapted = subset(fitness, generation == "G40")
adapted = merge(adapted,mean.fitness.gain.pop[,c("pop", "pop_ordered", "pop_random")])


mean.fitness.gain.rec=mean.fitness.gain.rec[order(mean.fitness.gain.rec$env, rev(mean.fitness.gain.rec$rec), decreasing = F),]
mean.fitness.gain.rec$xticks = xticks





FITNESS = list(ancetral.fitness = ancestor, 
               fitness.gain = adapted,
               mean.ancestral.fitness=mean.ancestral.fitness2,
               mean.fitness.gain.pop=mean.fitness.gain.pop,
               mean.fitness.gain.rec=mean.fitness.gain.rec)

save(FITNESS, file = "adaptation_fitness/fitnessGain.Rdata")

