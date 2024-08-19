
library(sommer)
library(ggplot2)
library(data.table)

prunethreshold = 0.9

PATH = "/Users/tomparee/Desktop/rec1EE_github/"
setwd(PATH)
source("./rec1EE_utils.R")

load(file = "heritability/phenotype/fertility.Rdata")
fertility = fertility[order(fertility$line),]


# Choose dataset = "CeMEEv2" or "markerSet1"
dataset = "CeMEEv2"
#dataset = "markerSet1"

#Fuctions: 

numtorom = function(vector_chromid){
  vector_chromid = as.character(as.roman(vector_chromid))
  vector_chromid[vector_chromid=="VI"]="X"
  vector_chromid = factor(vector_chromid, levels = c("I","II","III","IV","V","X"))
  return(vector_chromid)
}

romtonum = function(vector_chromid){
  vector_chromid[vector_chromid=="I"]=1
  vector_chromid[vector_chromid=="II"]=2
  vector_chromid[vector_chromid=="III"]=3
  vector_chromid[vector_chromid=="IV"]=4
  vector_chromid[vector_chromid=="V"]=5
  vector_chromid[vector_chromid=="X"]=6
  return(vector_chromid)
}

# standardized SNPs (all markers explain equal variance after standardization == low frequency markers explain more)  
# LD? every SNP (even those perfectly correlated)
# read speed and balding (https://www.ncbi.nlm.nih.gov/pubmed/25404112)
hgsm <- function(gmat){
  m = nrow(gmat) # markers
  gmat = scale(t(gmat))
  gmat[is.na(gmat)] <- 0
  # u = PCs, v = marker loadings, d = singular values for each PC
  # so that xx = svdx$u %*% diag(svdx$d) %*% t(svdx$v)
  # from Patterson 2006, Hoffman 2013, K = XX'= US(US)'
  X = tcrossprod(gmat) / m
  svdx <- svd(X)
  D=diag(sqrt(svdx$d))
  US = svdx$u %*% D
  K = tcrossprod(US)
  K <- K/mean(diag(K))
  rownames(K) <- colnames(K) <- rownames(gmat)
  K
}


doPrune = function(snps, X, r2=0.99, np=1, ...) {
  
  LDprune = function(snps, x, maxr2 = 0.99, window=5000, step=2000, randomSeed=F, maximise=NULL){
    # x=PxN matrix (single chromosome)
    # if randomSeed, choose tag marker within the window at random (rather than first)
    # if maximise, choose based on maximum value of a numeric vector
    P = nrow(x)
    X = t(x)
    colnames(X) = 1:P
    if(window > P) {
      wins = c(1, P+1)
    } else {
      wins = seq(1, P, window)
      wins = c(sort(c(wins, seq(step, P, window))), P+1)
    }
    keepers = outs = NULL
    for(i in 1:(length(wins)-1)){
      ix = wins[i]:(wins[i+1]-1)
      rx = cor(X[,ix])^2
      if(randomSeed) ix = sample(ix)
      if(!is.null(maximise)) {o=order(maximise[ix], decreasing = T); rx = rx[o,o]; ix=ix[o]}
      rx[upper.tri(rx, diag = T)] <- 0
      for(j in 1:nrow(rx)){
        if(!j %in% outs){
          jj = ix[which(rx[,j]>maxr2)]
          outs = append(outs, jj)
        }
      }
      outs = unique(outs)
      keepers = sort(unique(append(keepers, ix[!ix %in% outs])))
      if(!is.null(maximise)) if(max(maximise[keepers]) < max(maximise[outs])) cat(sprintf('error maximising: chrom %s window %s (%s-%s)\n', snps$chrom[1], i, snps$pos[min(ix)], snps$pos[max(ix)]))
    }
    list(snps[keepers,], x[keepers,])
  }
  
  while(1){
    n = nrow(snps)
    print(table(snps$chrom))
    lds <- parallel::mclapply(split(cbind(snps, X), snps$chrom), mc.cores = np, function(i) LDprune(i[,1:2], i[,-(1:2)], maxr2 = r2, ...))
    snps = do.call(rbind, lapply(lds,  '[[', 1))
    X = do.call(rbind, lapply(lds,  '[[', 2))
    if(nrow(snps)==n) {break; print(table(snps$chrom)); rownames(X) = rownames(snps) = NULL}
  }
  return(list(snps, X))
}

afsToMafs <- function(x) sapply(x, function(y) ifelse(y <= 0.5, y, 1-y))
filterMAF <- function(df, MAFgt=0, returnMAFs=F, returnSplit=F, header=4, lines=NULL) {
  # assumes chrom, pos, ref, alt in first four cols (header), all other cols genotypes (0/1)
  # if minMAF supplied, filters to > minMAF
  df <- data.frame(df)
  snps = df[,1:header]
  if(!is.null(lines)) {
    X <- df[,names(df) %in% lines]
    cat(sprintf('dropping from %s to %s lines (%s supplied)\n', ncol(df)-header, ncol(X), len(lines)))
  } else {
    X <- df[,-(1:header)]
  }
  nrils = ncol(X)
  counts = apply(X, 1, sum, na.rm=T)
  nnas = apply(X, 1, function(x) sum(!is.na(x)))
  mafs <- afsToMafs(counts/nnas)
  X <- X[mafs>MAFgt,]
  snps <- snps[mafs>MAFgt,]
  cat(sprintf('dropping %s markers with MAF <= %s, %s segregating\n', sum(mafs<=MAFgt), MAFgt, nrow(X)))
  if (returnMAFs){
    return(list(X, mafs))
  }
  if (returnSplit){
    return(list(snps, X))
  } else {
    return(cbind(snps, X))
  }
}


###############################################
#### Import the data ##########################

if(dataset == "markerSet1"){
  
  #gt <- fread('~/Documents/cemee/genotypes/WS220_CeMEEv2_markerSet1.csv.gz')
  #gmap <- fread('~/Documents/cemee/genotypes/WS220.F2.geneticMap.txt.gz')
  # add genetic distance
  #snps = merge(gmap, gt[, 1:2], sort=F)
  
  ## X is our transposed NxP (lines x markers) matrix
  #X = t(gt[,-(1:4)])
  # remove fixed sites
  #afs = apply(X, 2, sum)/nrow(X)
  #summary(afs)
  #X = X[,afs>0 & afs<1]
  #snps = snps[afs>0 & afs<1,]
  
  ## fertility BLUPs (from a mixed model using all replicate data)
  #load('~/Documents/cemee/phenotypes/fertility.rda', verbose = T)
  #table(ecoefs$env)
  
  #save(snps, gt, ecoefs, file = '~/Documents/cemee/rec1/geno_pheno_for_mapping.rda')
  
  
  load("geno_pheno_for_mapping.rda")
  gt = as.data.frame(gt)
  snps = gt[,1:4]
}


if(dataset == "CeMEEv2"){

  #Import genotype
  # (!) "CeMEEv2_RIL_geno.csv.gz" & CeMEEv2_RIL_snps_ws220.csv.gz need to be downloaded from:
  # https://github.com/lukemn/cemee/tree/gh-pages/cemee_v2_resource/supplement
  gt = as.data.frame(fread('/Users/tomparee/Documents/Documents - MacBook Pro de tom/cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_RIL_geno.csv.gz'))
  snps = as.data.frame(fread('/Users/tomparee/Documents/Documents - MacBook Pro de tom/cemee_LMN/GWAS/cemee_v2_resource/supplement/CeMEEv2_RIL_snps_ws220.csv.gz'))
  snps$chrom = tstrsplit(snps$chrom, "_")[[2]]
  snps$chrom = romtonum(snps$chrom)
  
  gt = filterMAF(cbind(snps, gt), MAFgt = 0.05)
  snps = gt[,1:4]
  gt = gt[,-(1:4)]
}

gt = gt[,colnames(gt) %in% fertility$line]
# snps$afm = apply(gt, 1, function(x) mean(x, na.rm=T))
# snps %>% ggplot(aes(pos, afm)) + geom_point(stroke=0, alpha=0.2) + facet_grid(.~chrom)

gt[gt==0.5] = NA # does not change much



afs = apply(gt, 1, function(x) sum(x, na.rm = T))/ncol(gt)
nas = apply(gt, 1, function(x) sum(is.na(x))==0)
ix = afs>0 & afs<1
# ix = afs>0 & afs<1 & nas
gt = gt[ix,]
snps = snps[ix,]

# mean impute missing data
gt = t(apply(gt, 1, function(x) {x[is.na(x)] = mean(x, na.rm=T); x}))

X = t(gt) ## X is our transposed NxP (lines x markers) matrix
dim(X)


gtld <- doPrune(snps[,1:2], gt, r2 = prunethreshold, np = 6)
snpsld = merge(gtld[[1]], snps, sort=F)
Xld = t(gtld[[2]])
dim(Xld)


snpsld$afm = apply(Xld, 2, function(x) mean(x, na.rm=T))
#ggplot(snpsld, aes(pos, afm)) + geom_point(stroke=0, alpha=0.2) + facet_grid(.~chrom)

# x = cor(Xld[,snpsld$chrom==1])^2
# image(x)
# y = x[upper.tri(x)]

#########################################################
##  add rrdiff (recombination difference; delta r) ####
load("recombinationMaps/recmaps.Rdata", verbose = T)
winsizeMb = 1
rrdiffthreshold = 0
snpsld = do.call(rbind, lapply(split(snpsld,snpsld$chrom), function(x){
  #x=split(snpsld,snpsld$chrom)[[3]]
  x = x[order(x$pos),]
  chr = x$chrom[1]
  
  ppos1 = x$pos - (1e6*winsizeMb)/2 
  ppos2 = x$pos + (1e6*winsizeMb)/2 
  
  b1 = ppos1 < 1
  b2 = ppos2 > max(x$pos)
  ppos2[b1] = ppos2[b1] - ppos1[b1]
  ppos1[b1] = 1
  
  ppos1[b2] = ppos1[b2] - (ppos2[b2] - max(x$pos))
  ppos2[b2] = max(x$pos)
  
  
  gwt = subset(recmaps, rec=='wt' & chrom == chr)
  gmut = subset(recmaps, rec=='mut' & chrom == chr)
  
  gwt1 = approx(x = gwt$pos, y = gwt$genetic, xout = ppos1)$y
  gwt2 = approx(x = gwt$pos, y = gwt$genetic, xout = ppos2)$y
  gmut1 = approx(x = gmut$pos, y = gmut$genetic, xout = ppos1)$y
  gmut2 = approx(x = gmut$pos, y = gmut$genetic, xout = ppos2)$y
  
  #x$log2rrchange = log( (gmut2-gmut1)/(gwt2-gwt1), base=2)
  pdist = ppos2-ppos1
  x$rrdiff = ((gmut2-gmut1)/pdist)-((gwt2-gwt1)/pdist)
  
  x$rec1effect = sign(x$rrdiff)
  x$rec1effect[x$rrdiff < rrdiffthreshold & x$rrdiff > -rrdiffthreshold] = 0
  x
}))


##########################################
### Calculate h2 #########################


lines = rownames(Xld)
phe = fertility[match(lines,fertility$line),]
phe$id1 = phe$id2 = phe$line

phe = subset(phe, line %in% lines)
phe = phe[order(phe$line),]
Xld = Xld[lines %in% phe$line,]
Xld = Xld[order(rownames(Xld)),]

# define the two regions: more or less recombination in mutant 
more = which(snpsld$rec1effect==1)
less = which(snpsld$rec1effect==-1)

hr1 <- hgsm(t(Xld[,more]))
hr2 <- hgsm(t(Xld[,less]))
#data.frame(x=as.numeric(hr1[upper.tri(hr1)])) %>% ggplot() + geom_histogram(aes(x))
#data.frame(x=as.numeric(hr2[upper.tri(hr2)])) %>% ggplot() + geom_histogram(aes(x))

# check against alternative additive GSM
#hr1.2 = A.mat(Xld[,more])
#hr2.2 <- A.mat(Xld[,less])
##data.frame(a=as.numeric(hr1), b= as.numeric(hr1.2)) %>% 
##  ggplot() + geom_point(aes(a, b))


# fit the model
require(sommer)
fitld <- mmer(fertility~1, random = ~vs(id1, Gu=hr1) + vs(id2, Gu=hr2), rcov = ~units, data=phe, date.warning = T, verbose = T)
#fitld <- mmer(fertility~1, random = ~vs(id1, Gu=hr1.2) + vs(id2, Gu=hr2.2), rcov = ~units, data=phe, date.warning = T, verbose = T)

vcs = unlist(fitld$sigma_scaled); vcs/sum(vcs)
h2 = data.frame(region = c("increased.rr", "reduced.rr", "all"),
                  h2=vcs/sum(vcs), 
                  SE =sqrt(diag(unlist(fitld$sigmaSE)))/sum(vcs),
                  nSNP = c(length(more), length(less), length(more)+length(less)))

save(h2,file = paste0("heritability/h2~deltaR_", prunethreshold, ".Rdata"))


#region        h2         SE  nSNP
#u:id1 increased.rr 0.1268435 0.01434275 11404
#u:id2   reduced.rr 0.2394717 0.01934435 36590
#units          all 0.6336848 0.02251489 47994