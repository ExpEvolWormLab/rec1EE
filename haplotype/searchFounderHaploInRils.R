PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

includeEEV1403AsFounder = F # if T, EEV1403 is considered as a founder 
# Note EEV1403 and EEV1401 have the same genetic background (only differ by cripr in rec-1)
# here it is referred as EEV1401

#####################################
#### UPLOAD FUNCTION AND DATA  ######

# Function:
source("haplotype/function_findhaplo_in_RILS.R")

# (!) you need to upload the different CeMEEv2 dataset from:
# https://github.com/lukemn/cemee/tree/gh-pages/cemee_v2_resource/supplement

foundersnps <- as.data.frame(fread('/Users/tomparee/Documents/Documents - MacBook Pro de tom/cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_founder_snps_ws245.csv.gz'))
founders <- as.data.frame(fread('/Users/tomparee/Documents/Documents - MacBook Pro de tom/cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_founder_geno.csv.gz'))


#### Merge the founders with highly similar haplotypes

# pairs = expand.grid(1:ncol(founders),1:ncol(founders)) # possible pairs of founders
# pairs = pairs[pairs$Var1<pairs$Var2,] #keep unique comb
# colnames(pairs) = c("founder1", "founder2")
# gdis = do.call(rbind, lapply(1:nrow(pairs), function(i){
#   
#  founder1 = founders[,pairs[i,1]]
#  founder2 = founders[,pairs[i,2]]
#  genotypedistance = sum(abs(founder1-founder2), na.rm=T)
#  
#  cbind(pairs[i,], gdis = genotypedistance, relgdis=genotypedistance/nrow(founders))
#   
# }))
# 
# gdis = gdis[order(gdis$gdis),]

rils = as.data.frame(fread("genotype/chrI_SMR_inbredlines_genotype_CeMEEv2_ws245.csv.gz"))
snps = as.data.frame(fread("genotype/chrI_SMR_inbredlines_snps_CeMEEv2_ws245.csv.gz"))


founders = founders[match(paste0(snps$chrom, snps$pos), paste0(foundersnps$chrom, foundersnps$pos)),]
#cemee = cemee[match(paste0(snps$chrom, snps$pos), paste0(cemeesnps$chrom, cemeesnps$pos)),]



if(includeEEV1403AsFounder==T){
  founders$EEV1403 = rils[,'EEV1401']
  rils = rils[,-which(colnames(rils)=='EEV1401')]
}




#######################################
############## RUN ####################
#######################################

ix = 1:ncol(rils)
haplo = do.call(rbind,parallel::mclapply(ix,mc.cores = 2, function(thisril){
  print(thisril)
  ril = rils[,thisril]
  ril[ril!=1 & ril!=0]=NA
  output = try(haplosearch(ril, founders=founders, snps))
  if(class(output)!='try-error'){
    output$rilname = colnames(rils)[thisril]
    return(output)
  }else{
    return(NULL)
  }
}))

#unique(haplo$rilname)

#save(haplo,file="haplotype/RILs_haplotypes_founders&EEV1403.Rdata")
#save(haplo,file="haplotype/RILs_haplotypes_founderOnly.Rdata")

#load(file="haplotype/RILs_haplotypes_foundersOnly.Rdata")
#load(file="haplotype/RILs_haplotypes_founders&EEV1403.Rdata")


### Put an arbitrary breakpoint where the haplotype overlaps

foundernames = colnames(haplo)[3:(ncol(haplo)-1)]

haplo$pos1 = snps$pos[haplo$whichsnp1]
haplo$pos2 = snps$pos[haplo$whichsnp2]

haplo = do.call(rbind,lapply(split(haplo, haplo$rilname), function(HAPLO){
  
  #HAPLO=split(haplo, haplo$rilname)[[5]]
  
  print(HAPLO$rilname[1])
  
  HAPLO = FuseHapSeparatedByOneSNP(HAPLO=HAPLO, posinfo = c("whichsnp", "pos"),
                                   foundercolnames=foundernames)
  
  
  if(nrow(HAPLO)>1){HAPLO = putBreakAtMidDistance(HAPLO=HAPLO,colnames.distance = c("pos1", "pos2"),round=T)}
  
  HAPLO$whichsnp1 = ceiling(approx(x=snps$pos, y=1:nrow(snps), xout=HAPLO$pos1)$y)
  HAPLO$whichsnp2 = floor(approx(x=snps$pos, y=1:nrow(snps), xout=HAPLO$pos2)$y)
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  overlap1snp = which(HAPLO$whichsnp1[2:nrow(HAPLO)] == HAPLO$whichsnp2[1:(nrow(HAPLO)-1)])+1
  HAPLO$whichsnp1[overlap1snp] = HAPLO$whichsnp1[overlap1snp]+1
  HAPLO[HAPLO$whichsnp1>HAPLO$whichsnp2,]
  HAPLO
}))

### Choose the most frequent founder when several founders are possible
haplo = KeepMostFrequentHaplo(haplotype=haplo, foundernames=foundernames)




#save(haplo,file="haplotype/RILs_haplotypes_foundersOnly_resolved.Rdata")
#save(haplo,file="haplotype/RILs_haplotypes_founders&EEV1403_resolved.Rdata")


