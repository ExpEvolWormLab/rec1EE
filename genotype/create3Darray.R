
PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

library(data.table)


################################
### Transform in a 3D array  ###


snpsname = "Rhom_poolseq_snps_CeMEEv2_ws245.csv.gz"
ADname = "Rhom_poolseq_AllelicDepth_CeMEEv2_ws245.csv.gz"

AD=NULL
snps=NULL
for(chr in c("I", "II", "III", "IV", 'V', 'X')){
  
   snpschr = fread(paste0("genotype/chr",chr,'_', snpsname))
   ADchr = fread(paste0("genotype/chr",chr,'_', ADname))
   
   snps = rbind(snps, snpschr)
   AD = rbind(AD, ADchr)
   
}

AD=as.data.frame(AD)
snps = as.data.frame(snps)

# split AD in two df, one for ref count, and the second for alt count
ref = do.call(cbind, lapply(1:ncol(AD), function(i){ as.numeric(tstrsplit(AD[,i],':')[[1]]) }))
alt = do.call(cbind, lapply(1:ncol(AD), function(i){ as.numeric(tstrsplit(AD[,i],':')[[2]]) }))
colnames(ref)=colnames(alt)=colnames(AD)

#transfor alt & ref in 3D array

gen = unique( tstrsplit(colnames(ref), "_")[[2]])
gen = gen[order(as.numeric(tstrsplit(gen, "G")[[2]]))]

pop =  unique( tstrsplit(colnames(ref), "_")[[1]])
pop = pop[order(pop)]

marker = paste0(snps$chrom, '_' ,snps$pos)

ref = lapply(gen, function(x){
  x = ref[,grepl(x, colnames(ref))]
  x = x[,match(pop, tstrsplit(colnames(x), '_')[[1]] )]
  colnames(x) = pop
  x
})

names(ref) = gen
ref = array(unlist(ref), dim = c(length(marker), length(pop), length(gen)), dimnames = list(marker, pop, gen))


alt = lapply(gen, function(x){
  x = alt[,grepl(x, colnames(alt))]
  x = x[,match(pop, tstrsplit(colnames(x), '_')[[1]] )]
  colnames(x) = pop
  x
})

names(alt) = gen

alt = array(unlist(alt), dim = c(length(marker), length(pop), length(gen)), dimnames = list(marker, pop, gen))

# Save the data (unfiltered)
save(snps, alt, ref, file = paste0("genotype/", "Rhom_poolsed_unfiltered.Rdata"))



#Filter data (rm outlier population at concerned timepoints + filterout low depth SNV)

outlier = rbind(data.frame(pop = "LR1", gen = c("G30", "G33", "G41")),
                data.frame(pop = "SHR4", gen = c("G24","G30", "G33", "G41")))

for(i in 1:nrow(outlier)){
  ref[, outlier[i,"pop"], outlier[i,"gen"]]=NA
  alt[, outlier[i,"pop"], outlier[i,"gen"]]=NA
}


depth = ref + alt
depth = apply(depth, 1, sum, na.rm=T)
hist(depth, breaks = 100)

lowdepth = depth<quantile(depth,0.05)
fixed = apply(alt, 1, sum, na.rm=T) == 0 | apply(ref, 1, sum, na.rm=T) ==0



ref = ref[!lowdepth & !fixed,,]
alt = alt[!lowdepth & !fixed,,]
snps = snps[!lowdepth & !fixed,]

# Save the data (filtered)
save(snps, alt, ref, file = paste0("genotype/", "Rhom_poolseq.Rdata"))


